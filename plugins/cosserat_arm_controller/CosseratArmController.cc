#include "CosseratArmController.hh"

#include <algorithm>
#include <chrono>
#include <cmath>

#include <gz/plugin/Register.hh>
#include <gz/common/Console.hh>

#include <gz/sim/Joint.hh>
#include <gz/sim/Link.hh>

using namespace cosserat_arm_controller;

namespace
{
/// Extracts the 2-DOF (theta_x, theta_y) universal-joint angles from a
/// relative rotation _dR = Rx(theta_x) * Ry(theta_y) - the same
/// intrinsic-X-then-Y composition SDF/DART use for a <joint type="universal">
/// with axis=(1,0,0), axis2=(0,1,0). Derived from expanding
/// Rx(tx)*Ry(ty) symbolically:
///   dR(1,1) = cos(tx),           dR(2,1) = sin(tx)
///   dR(0,0) = cos(ty),           dR(0,2) = sin(ty)
std::pair<double, double> ExtractXYEuler(const cosserat_arm_controller::Mat3 &_dR)
{
  const double thetaX = std::atan2(_dR.m[2][1], _dR.m[1][1]);
  const double thetaY = std::atan2(_dR.m[0][2], _dR.m[0][0]);
  return {thetaX, thetaY};
}
}  // namespace

//////////////////////////////////////////////////
void CosseratArmController::Configure(
    const gz::sim::Entity &_entity,
    const std::shared_ptr<const sdf::Element> &_sdf,
    gz::sim::EntityComponentManager &_ecm,
    gz::sim::EventManager & /*_eventMgr*/)
{
  this->model = gz::sim::Model(_entity);
  if (!this->model.Valid(_ecm))
  {
    gzerr << "CosseratArmController must be attached to a model entity."
          << std::endl;
    return;
  }

  if (auto elem = _sdf->FindElement("joint_name"))
  {
    while (elem)
    {
      this->jointNames.push_back(elem->Get<std::string>());
      elem = elem->GetNextElement("joint_name");
    }
  }
  for (const auto &name : this->jointNames)
  {
    gz::sim::Entity joint = this->model.JointByName(_ecm, name);
    if (joint == gz::sim::kNullEntity)
    {
      gzerr << "CosseratArmController: could not find joint [" << name << "]"
            << std::endl;
    }
    else
    {
      gz::sim::Joint(joint).EnablePositionCheck(_ecm, true);
      gz::sim::Joint(joint).EnableVelocityCheck(_ecm, true);
    }
    this->jointEntities.push_back(joint);
  }
  this->jointTargets.assign(this->jointEntities.size(), {0.0, 0.0});

  CosseratRodSolver::Params params;
  if (_sdf->HasElement("rod_length"))
    params.length = _sdf->Get<double>("rod_length");
  if (_sdf->HasElement("rod_outer_radius"))
    params.outerRadius = _sdf->Get<double>("rod_outer_radius");
  if (_sdf->HasElement("rod_inner_radius"))
    params.innerRadius = _sdf->Get<double>("rod_inner_radius");
  if (_sdf->HasElement("youngs_modulus"))
    params.youngsModulus = _sdf->Get<double>("youngs_modulus");
  if (_sdf->HasElement("shear_modulus"))
    params.shearModulus = _sdf->Get<double>("shear_modulus");
  if (_sdf->HasElement("solidity_fill_factor"))
    params.solidityFillFactor = _sdf->Get<double>("solidity_fill_factor");
  if (_sdf->HasElement("rod_density"))
    params.density = _sdf->Get<double>("rod_density");
  if (_sdf->HasElement("pitch_radius"))
    params.tendonPitchRadius = _sdf->Get<double>("pitch_radius");
  if (_sdf->HasElement("distribute_tendon_load_continuously"))
    params.distributeTendonLoadContinuously =
        _sdf->Get<bool>("distribute_tendon_load_continuously");
  if (_sdf->HasElement("use_guided_tendons"))
    this->useGuidedTendons = _sdf->Get<bool>("use_guided_tendons");
  if (this->useGuidedTendons)
  {
    // SolveGuided() needs actual discrete notches to compute chord
    // directions between (see CosseratRodSolver.hpp's SolveGuided() and
    // Params::tendonNotchFractions comments) - thirds, matching the
    // physical hardware's 3 internal routing notches.
    params.tendonNotchFractions = {1.0 / 3.0, 2.0 / 3.0, 1.0};
    params.distributeTendonLoadContinuously = false;
  }
  this->restLength = params.length;
  this->solver = std::make_unique<CosseratRodSolver>(params);

  if (_sdf->HasElement("tendon_stiffness"))
    this->tendonStiffness = _sdf->Get<double>("tendon_stiffness");
  if (_sdf->HasElement("max_tension"))
    this->maxTension = _sdf->Get<double>("max_tension");
  if (_sdf->HasElement("tendon_baseline"))
    this->tendonBaseline = _sdf->Get<double>("tendon_baseline");
  if (_sdf->HasElement("position_p_gain"))
    this->positionPGain = _sdf->Get<double>("position_p_gain");
  if (_sdf->HasElement("position_d_gain"))
    this->positionDGain = _sdf->Get<double>("position_d_gain");
  if (_sdf->HasElement("effort_limit"))
    this->effortLimit = _sdf->Get<double>("effort_limit");
  if (_sdf->HasElement("tendon_topic"))
    this->tendonTopic = _sdf->Get<std::string>("tendon_topic");
  if (_sdf->HasElement("tip_pose_topic"))
    this->tipPoseTopic = _sdf->Get<std::string>("tip_pose_topic");
  if (_sdf->HasElement("command_filter_omega"))
    this->commandFilterOmega = _sdf->Get<double>("command_filter_omega");
  if (_sdf->HasElement("max_tendon_rate_m_s"))
    this->maxTendonRateMPerS = _sdf->Get<double>("max_tendon_rate_m_s");
  if (_sdf->HasElement("solve_rate_hz"))
    this->solveRateHz = _sdf->Get<double>("solve_rate_hz");

  this->shaperL1.SetOmega(this->commandFilterOmega);
  this->shaperL2.SetOmega(this->commandFilterOmega);
  this->shaperL3.SetOmega(this->commandFilterOmega);
  this->shaperL1.SetMaxRate(this->maxTendonRateMPerS);
  this->shaperL2.SetMaxRate(this->maxTendonRateMPerS);
  this->shaperL3.SetMaxRate(this->maxTendonRateMPerS);
  this->shaperL1.Reset(this->restLength);
  this->shaperL2.Reset(this->restLength);
  this->shaperL3.Reset(this->restLength);

  if (_sdf->HasElement("tip_link"))
  {
    const std::string tipLinkName = _sdf->Get<std::string>("tip_link");
    this->tipLinkEntity = this->model.LinkByName(_ecm, tipLinkName);
    if (this->tipLinkEntity == gz::sim::kNullEntity)
    {
      gzerr << "CosseratArmController: could not find tip link ["
            << tipLinkName << "]" << std::endl;
    }
  }

  this->tipPosePub = this->node.Advertise<gz::msgs::Pose>(this->tipPoseTopic);
  this->computeTimePub =
      this->node.Advertise<gz::msgs::Double>(this->tendonTopic + "/compute_time_ms");
  this->node.Subscribe(this->tendonTopic, &CosseratArmController::OnTendonCmd,
                        this);

  gzmsg << "CosseratArmController configured with "
        << this->jointEntities.size() << " joints, rod length "
        << params.length << "m: command topic [" << this->tendonTopic
        << "], tip pose topic [" << this->tipPoseTopic << "]" << std::endl;
}

//////////////////////////////////////////////////
void CosseratArmController::OnTendonCmd(const gz::msgs::Vector3d &_msg)
{
  // Only updates the CommandShapers' targets - PreUpdate() advances them
  // every physics tick and RunSolve() (throttled to solveRateHz) consumes
  // the FILTERED lengths, not this raw command directly. See class header
  // comment.
  std::lock_guard<std::mutex> lock(this->mutex);
  this->shaperL1.SetTarget(_msg.x());
  this->shaperL2.SetTarget(_msg.y());
  this->shaperL3.SetTarget(_msg.z());
}

//////////////////////////////////////////////////
void CosseratArmController::RunSolve(double _l1, double _l2, double _l3)
{
  // Every tendon carries at least tendonBaseline_ N even when un-shortened
  // (l_i >= restLength) - a real cable-driven arm keeps its tendons taut
  // rather than fully slack, and (empirically, see tendonBaseline_'s
  // comment) this also moves the solver's operating point onto a much
  // larger clean/reversal-free branch of its tension-vs-bend response.
  std::vector<double> tensions(3);
  tensions[0] = std::clamp(
      this->tendonBaseline +
          this->tendonStiffness * std::max(0.0, this->restLength - _l1),
      0.0, this->maxTension);
  tensions[1] = std::clamp(
      this->tendonBaseline +
          this->tendonStiffness * std::max(0.0, this->restLength - _l2),
      0.0, this->maxTension);
  tensions[2] = std::clamp(
      this->tendonBaseline +
          this->tendonStiffness * std::max(0.0, this->restLength - _l3),
      0.0, this->maxTension);

  const bool converged = this->useGuidedTendons
      ? this->solver->SolveGuided(tensions)
      : this->solver->Solve(tensions);
  if (!converged)
  {
    gzwarn << "CosseratArmController: BVP shooting solve did not converge "
              "for filtered tendon lengths (" << _l1 << ", " << _l2 << ", "
           << _l3 << ") - using best-effort iterate." << std::endl;
  }

  const int n = static_cast<int>(this->jointEntities.size());
  auto samples = this->solver->SampleSegmentCenters(n);

  std::vector<std::pair<double, double>> targets(n, {0.0, 0.0});
  cosserat_arm_controller::Mat3 prevR = cosserat_arm_controller::Mat3::Identity();
  for (int i = 0; i < n && i < static_cast<int>(samples.size()); ++i)
  {
    const cosserat_arm_controller::Mat3 dR = prevR.Transpose() * samples[i].R;
    auto [thetaX, thetaY] = ExtractXYEuler(dR);
    // The solver is a standalone numerical component with its own abstract
    // local frame (R(0)=Identity at the base, +z growth direction); it has
    // no notion of the model's actual world orientation. The old
    // {thetaX, -thetaY} correction here (a single Y-negation) was live-
    // tested (single-tendon l1/l2/l3 commands on both continuum_arm_pcc and
    // continuum_arm_cosserat side by side, see side_by_side_demo.py) and
    // found to leave the two models bending toward azimuths offset by
    // ~180deg from TendonArmController::UpdateTargets()'s convention (see
    // TendonArmController.cc) -- i.e. that mapping was the exact negation
    // of the correct one. Dropping the negation, a plain swap
    // {thetaY, thetaX}, was live-verified to bring both models' measured
    // tip azimuth within ~0.2deg of each other for all three tendons.
    targets[i] = {thetaY, thetaX};
    prevR = samples[i].R;
  }

  std::lock_guard<std::mutex> lock(this->mutex);
  this->jointTargets = targets;
}

//////////////////////////////////////////////////
void CosseratArmController::PreUpdate(
    const gz::sim::UpdateInfo &_info, gz::sim::EntityComponentManager &_ecm)
{
  if (_info.paused)
    return;

  const double dt = std::chrono::duration<double>(_info.dt).count();

  double filteredL1, filteredL2, filteredL3;
  bool doSolve = false;
  {
    std::lock_guard<std::mutex> lock(this->mutex);
    this->shaperL1.Update(dt);
    this->shaperL2.Update(dt);
    this->shaperL3.Update(dt);
    filteredL1 = this->shaperL1.Value();
    filteredL2 = this->shaperL2.Value();
    filteredL3 = this->shaperL3.Value();

    this->solveAccumulator += dt;
    const double solvePeriod = 1.0 / std::max(this->solveRateHz, 1.0);
    if (this->solveAccumulator >= solvePeriod)
    {
      this->solveAccumulator = std::fmod(this->solveAccumulator, solvePeriod);
      // Only actually re-solve if the filtered input has moved since the
      // last solve - see solveDeadBandM's declaration comment for why
      // continuously re-solving a static input is harmful, not just
      // wasteful.
      const double moved = std::max({std::abs(filteredL1 - this->lastSolvedL1),
                                      std::abs(filteredL2 - this->lastSolvedL2),
                                      std::abs(filteredL3 - this->lastSolvedL3)});
      if (!this->haveSolved || moved > this->solveDeadBandM)
      {
        doSolve = true;
        this->lastSolvedL1 = filteredL1;
        this->lastSolvedL2 = filteredL2;
        this->lastSolvedL3 = filteredL3;
        this->haveSolved = true;
      }
    }
  }

  if (doSolve)
  {
    const auto solveStart = std::chrono::steady_clock::now();
    this->RunSolve(filteredL1, filteredL2, filteredL3);
    const double solveMs = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - solveStart).count();
    gz::msgs::Double computeMsg;
    computeMsg.set_data(solveMs);
    this->computeTimePub.Publish(computeMsg);
  }

  std::vector<std::pair<double, double>> targets;
  {
    std::lock_guard<std::mutex> lock(this->mutex);
    targets = this->jointTargets;
  }

  for (std::size_t i = 0; i < this->jointEntities.size(); ++i)
  {
    const gz::sim::Entity jointEntity = this->jointEntities[i];
    if (jointEntity == gz::sim::kNullEntity || i >= targets.size())
      continue;

    gz::sim::Joint joint(jointEntity);
    auto position = joint.Position(_ecm);
    auto velocity = joint.Velocity(_ecm);
    if (!position.has_value() || position->size() < 2 ||
        !velocity.has_value() || velocity->size() < 2)
    {
      continue;
    }

    const double targetX = targets[i].first;
    const double targetY = targets[i].second;

    double errorX = targetX - (*position)[0];
    double errorY = targetY - (*position)[1];
    double effortX = this->positionPGain * errorX -
                      this->positionDGain * (*velocity)[0];
    double effortY = this->positionPGain * errorY -
                      this->positionDGain * (*velocity)[1];
    effortX = std::clamp(effortX, -this->effortLimit, this->effortLimit);
    effortY = std::clamp(effortY, -this->effortLimit, this->effortLimit);
    joint.SetForce(_ecm, {effortX, effortY});
  }
}

//////////////////////////////////////////////////
void CosseratArmController::PostUpdate(
    const gz::sim::UpdateInfo & /*_info*/,
    const gz::sim::EntityComponentManager &_ecm)
{
  if (this->tipLinkEntity == gz::sim::kNullEntity)
    return;

  auto tipPose = gz::sim::Link(this->tipLinkEntity).WorldPose(_ecm);
  if (!tipPose.has_value())
    return;

  gz::msgs::Pose msg;
  msg.mutable_position()->set_x(tipPose->Pos().X());
  msg.mutable_position()->set_y(tipPose->Pos().Y());
  msg.mutable_position()->set_z(tipPose->Pos().Z());
  msg.mutable_orientation()->set_w(tipPose->Rot().W());
  msg.mutable_orientation()->set_x(tipPose->Rot().X());
  msg.mutable_orientation()->set_y(tipPose->Rot().Y());
  msg.mutable_orientation()->set_z(tipPose->Rot().Z());
  this->tipPosePub.Publish(msg);
}

GZ_ADD_PLUGIN(CosseratArmController,
              gz::sim::System,
              gz::sim::ISystemConfigure,
              gz::sim::ISystemPreUpdate,
              gz::sim::ISystemPostUpdate)
