#include "RbsnArmController.hh"

#include <algorithm>
#include <chrono>
#include <cmath>

#include <gz/plugin/Register.hh>
#include <gz/common/Console.hh>

#include <gz/sim/Joint.hh>

using namespace rbsn_arm_controller;

namespace
{
// sqrt(3)/2, used by the tendon-routing Jacobian below (3 tendons at
// 0/120/240 deg on the pitch circle).
constexpr double kSqrt3Over2 = 0.8660254037844386;
}  // namespace

//////////////////////////////////////////////////
void RbsnArmController::Configure(
    const gz::sim::Entity &_entity,
    const std::shared_ptr<const sdf::Element> &_sdf,
    gz::sim::EntityComponentManager &_ecm,
    gz::sim::EventManager & /*_eventMgr*/)
{
  this->model = gz::sim::Model(_entity);
  if (!this->model.Valid(_ecm))
  {
    gzerr << "RbsnArmController must be attached to a model entity."
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
      gzerr << "RbsnArmController: could not find joint [" << name << "]"
            << std::endl;
    }
    else
    {
      gz::sim::Joint(joint).EnablePositionCheck(_ecm, true);
      gz::sim::Joint(joint).EnableVelocityCheck(_ecm, true);
    }
    this->jointEntities.push_back(joint);
  }

  if (_sdf->HasElement("pitch_radius"))
    this->pitchRadius = _sdf->Get<double>("pitch_radius");
  if (_sdf->HasElement("rest_length"))
    this->restLength = _sdf->Get<double>("rest_length");
  if (_sdf->HasElement("tendon_stiffness"))
    this->tendonStiffness = _sdf->Get<double>("tendon_stiffness");
  if (_sdf->HasElement("max_tension"))
    this->maxTension = _sdf->Get<double>("max_tension");
  if (_sdf->HasElement("tendon_baseline"))
    this->tendonBaseline = _sdf->Get<double>("tendon_baseline");
  if (_sdf->HasElement("k1"))
    this->k1 = _sdf->Get<double>("k1");
  if (_sdf->HasElement("k2"))
    this->k2 = _sdf->Get<double>("k2");
  if (_sdf->HasElement("c_theta"))
    this->cTheta = _sdf->Get<double>("c_theta");
  if (_sdf->HasElement("effort_limit"))
    this->effortLimit = _sdf->Get<double>("effort_limit");
  if (_sdf->HasElement("tendon_topic"))
    this->tendonTopic = _sdf->Get<std::string>("tendon_topic");
  if (_sdf->HasElement("command_filter_omega"))
    this->commandFilterOmega = _sdf->Get<double>("command_filter_omega");
  if (_sdf->HasElement("max_tendon_rate_m_s"))
    this->maxTendonRateMPerS = _sdf->Get<double>("max_tendon_rate_m_s");

  this->shaperL1.SetOmega(this->commandFilterOmega);
  this->shaperL2.SetOmega(this->commandFilterOmega);
  this->shaperL3.SetOmega(this->commandFilterOmega);
  this->shaperL1.SetMaxRate(this->maxTendonRateMPerS);
  this->shaperL2.SetMaxRate(this->maxTendonRateMPerS);
  this->shaperL3.SetMaxRate(this->maxTendonRateMPerS);
  this->shaperL1.Reset(this->restLength);
  this->shaperL2.Reset(this->restLength);
  this->shaperL3.Reset(this->restLength);

  this->computeTimePub =
      this->node.Advertise<gz::msgs::Double>(this->tendonTopic + "/compute_time_ms");
  this->node.Subscribe(this->tendonTopic, &RbsnArmController::OnTendonCmd,
                        this);

  gzmsg << "RbsnArmController configured with " << this->jointEntities.size()
        << " joints: command topic [" << this->tendonTopic << "]"
        << std::endl;
}

//////////////////////////////////////////////////
void RbsnArmController::OnTendonCmd(const gz::msgs::Vector3d &_msg)
{
  std::lock_guard<std::mutex> lock(this->mutex);
  this->shaperL1.SetTarget(_msg.x());
  this->shaperL2.SetTarget(_msg.y());
  this->shaperL3.SetTarget(_msg.z());
}

//////////////////////////////////////////////////
void RbsnArmController::PreUpdate(
    const gz::sim::UpdateInfo &_info, gz::sim::EntityComponentManager &_ecm)
{
  if (_info.paused)
    return;

  const auto computeStart = std::chrono::steady_clock::now();

  const double dt = std::chrono::duration<double>(_info.dt).count();

  double l1, l2, l3;
  {
    std::lock_guard<std::mutex> lock(this->mutex);
    this->shaperL1.Update(dt);
    this->shaperL2.Update(dt);
    this->shaperL3.Update(dt);
    l1 = this->shaperL1.Value();
    l2 = this->shaperL2.Value();
    l3 = this->shaperL3.Value();
  }

  // Tendon length -> pulling tension (a tendon can only pull, never push),
  // same spring-clutch model continuum_arm_cosserat uses.
  const double t1 = std::clamp(
      this->tendonBaseline + this->tendonStiffness * std::max(0.0, this->restLength - l1),
      0.0, this->maxTension);
  const double t2 = std::clamp(
      this->tendonBaseline + this->tendonStiffness * std::max(0.0, this->restLength - l2),
      0.0, this->maxTension);
  const double t3 = std::clamp(
      this->tendonBaseline + this->tendonStiffness * std::max(0.0, this->restLength - l3),
      0.0, this->maxTension);

  // Tendon-routing Jacobian transpose: converts pulling tensions directly
  // into a generalized actuation TORQUE (tau_x, tau_y), the force-domain
  // dual of TendonArmController::UpdateTargets()'s length-to-angle
  // kinematic map. Derivation (tendons at pitch radius d, azimuths
  // 0/120/240 deg): the forward kinematic Jacobian relating small joint
  // rotations to tendon length changes is
  //   dl1 = -d*dThetaX
  //   dl2 = +0.5*d*dThetaX - (sqrt3/2)*d*dThetaY
  //   dl3 = +0.5*d*dThetaX + (sqrt3/2)*d*dThetaY
  // (this is exactly the algebraic inverse of UpdateTargets()'s
  // theta = -(2/(3d))*(...) map - substituting confirms round-trip
  // consistency). By the standard robot-Jacobian force/torque duality
  // (power balance: -T^T*d(dl)/dt = tau^T*d(theta)/dt for all joint
  // velocities), tau = -J^T*T, which works out to:
  const double totalTauX = this->pitchRadius * (t1 - 0.5 * t2 - 0.5 * t3);
  const double totalTauY = this->pitchRadius * kSqrt3Over2 * (t2 - t3);

  const std::size_t n = std::max<std::size_t>(this->jointEntities.size(), 1);
  const double actuationX = totalTauX / static_cast<double>(n);
  const double actuationY = totalTauY / static_cast<double>(n);

  for (const auto &jointEntity : this->jointEntities)
  {
    if (jointEntity == gz::sim::kNullEntity)
      continue;

    gz::sim::Joint joint(jointEntity);
    auto position = joint.Position(_ecm);
    auto velocity = joint.Velocity(_ecm);
    if (!position.has_value() || position->size() < 2 ||
        !velocity.has_value() || velocity->size() < 2)
    {
      continue;
    }

    const double thetaX = (*position)[0];
    const double thetaY = (*position)[1];
    const double thetaXDot = (*velocity)[0];
    const double thetaYDot = (*velocity)[1];

    const double springX = this->k1 * thetaX + this->k2 * thetaX * thetaX * thetaX;
    const double springY = this->k1 * thetaY + this->k2 * thetaY * thetaY * thetaY;
    const double dampX = this->cTheta * thetaXDot;
    const double dampY = this->cTheta * thetaYDot;

    double effortX = actuationX - springX - dampX;
    double effortY = actuationY - springY - dampY;
    effortX = std::clamp(effortX, -this->effortLimit, this->effortLimit);
    effortY = std::clamp(effortY, -this->effortLimit, this->effortLimit);
    joint.SetForce(_ecm, {effortX, effortY});
  }

  const double computeMs = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - computeStart).count();
  gz::msgs::Double computeMsg;
  computeMsg.set_data(computeMs);
  this->computeTimePub.Publish(computeMsg);
}

GZ_ADD_PLUGIN(RbsnArmController,
              gz::sim::System,
              gz::sim::ISystemConfigure,
              gz::sim::ISystemPreUpdate)
