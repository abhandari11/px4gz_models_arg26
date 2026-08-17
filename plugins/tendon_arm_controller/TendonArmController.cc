#include "TendonArmController.hh"

#include <algorithm>
#include <chrono>
#include <cmath>

#include <gz/plugin/Register.hh>
#include <gz/common/Console.hh>

#include <gz/sim/Joint.hh>

using namespace tendon_arm_controller;

//////////////////////////////////////////////////
void TendonArmController::Configure(
    const gz::sim::Entity &_entity,
    const std::shared_ptr<const sdf::Element> &_sdf,
    gz::sim::EntityComponentManager &_ecm,
    gz::sim::EventManager & /*_eventMgr*/)
{
  this->model = gz::sim::Model(_entity);
  if (!this->model.Valid(_ecm))
  {
    gzerr << "TendonArmController must be attached to a model entity."
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
      gzerr << "TendonArmController: could not find joint [" << name << "]"
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
  if (_sdf->HasElement("position_p_gain"))
    this->positionPGain = _sdf->Get<double>("position_p_gain");
  if (_sdf->HasElement("position_d_gain"))
    this->positionDGain = _sdf->Get<double>("position_d_gain");
  if (_sdf->HasElement("effort_limit"))
    this->effortLimit = _sdf->Get<double>("effort_limit");
  if (_sdf->HasElement("tendon_topic"))
    this->tendonTopic = _sdf->Get<std::string>("tendon_topic");
  if (_sdf->HasElement("angles_topic"))
    this->anglesTopic = _sdf->Get<std::string>("angles_topic");
  if (_sdf->HasElement("command_filter_omega"))
    this->commandFilterOmega = _sdf->Get<double>("command_filter_omega");
  if (_sdf->HasElement("max_tendon_rate_m_s"))
    this->maxTendonRateMPerS = _sdf->Get<double>("max_tendon_rate_m_s");

  // Zero tendon deflection (l1=l2=l3=restLength) until the first command.
  this->l1 = this->l2 = this->l3 = this->restLength;
  this->l1Cmd = this->l2Cmd = this->l3Cmd = this->restLength;
  this->shaperL1.SetOmega(this->commandFilterOmega);
  this->shaperL2.SetOmega(this->commandFilterOmega);
  this->shaperL3.SetOmega(this->commandFilterOmega);
  this->shaperL1.SetMaxRate(this->maxTendonRateMPerS);
  this->shaperL2.SetMaxRate(this->maxTendonRateMPerS);
  this->shaperL3.SetMaxRate(this->maxTendonRateMPerS);
  this->shaperL1.Reset(this->restLength);
  this->shaperL2.Reset(this->restLength);
  this->shaperL3.Reset(this->restLength);
  this->shapersInitialized = true;

  this->anglesPub =
      this->node.Advertise<gz::msgs::Vector3d>(this->anglesTopic);
  this->computeTimePub =
      this->node.Advertise<gz::msgs::Double>(this->tendonTopic + "/compute_time_ms");

  this->node.Subscribe(this->tendonTopic, &TendonArmController::OnTendonCmd,
                        this);

  gzmsg << "TendonArmController configured with " << this->jointEntities.size()
        << " joints: command topic [" << this->tendonTopic
        << "], angles topic [" << this->anglesTopic << "]" << std::endl;
}

//////////////////////////////////////////////////
void TendonArmController::UpdateTargets()
{
  const double d = this->pitchRadius;
  const double dl1 = this->l1 - this->restLength;
  const double dl2 = this->l2 - this->restLength;
  const double dl3 = this->l3 - this->restLength;
  const double kSqrt3Over2 = 0.8660254037844386;

  this->thetaX = -(2.0 / (3.0 * d)) * (dl1 - 0.5 * dl2 - 0.5 * dl3);
  this->thetaY =
      -(2.0 / (3.0 * d)) * (kSqrt3Over2 * dl2 - kSqrt3Over2 * dl3);

  const std::size_t n = std::max<std::size_t>(this->jointEntities.size(), 1);
  this->thetaXPerJoint = this->thetaX / static_cast<double>(n);
  this->thetaYPerJoint = this->thetaY / static_cast<double>(n);
}

//////////////////////////////////////////////////
void TendonArmController::OnTendonCmd(const gz::msgs::Vector3d &_msg)
{
  // Raw commands only set the CommandShapers' targets - PreUpdate() advances
  // them every physics tick and feeds the FILTERED length into
  // UpdateTargets(), giving the arm its smooth, zero-overshoot, ~2s servo
  // response instead of an instant target jump (see CommandShaper.hpp).
  std::lock_guard<std::mutex> lock(this->mutex);
  this->l1Cmd = _msg.x();
  this->l2Cmd = _msg.y();
  this->l3Cmd = _msg.z();
  this->shaperL1.SetTarget(this->l1Cmd);
  this->shaperL2.SetTarget(this->l2Cmd);
  this->shaperL3.SetTarget(this->l3Cmd);
  this->haveCommand = true;
}

//////////////////////////////////////////////////
void TendonArmController::PreUpdate(
    const gz::sim::UpdateInfo &_info, gz::sim::EntityComponentManager &_ecm)
{
  if (_info.paused)
    return;

  const auto computeStart = std::chrono::steady_clock::now();

  const double dt =
      std::chrono::duration<double>(_info.dt).count();

  double targetX, targetY;
  gz::msgs::Vector3d anglesMsg;
  {
    std::lock_guard<std::mutex> lock(this->mutex);
    this->shaperL1.Update(dt);
    this->shaperL2.Update(dt);
    this->shaperL3.Update(dt);
    this->l1 = this->shaperL1.Value();
    this->l2 = this->shaperL2.Value();
    this->l3 = this->shaperL3.Value();
    this->UpdateTargets();
    targetX = this->thetaXPerJoint;
    targetY = this->thetaYPerJoint;
    anglesMsg.set_x(this->thetaX);
    anglesMsg.set_y(this->thetaY);
    anglesMsg.set_z(0.0);
  }
  this->anglesPub.Publish(anglesMsg);

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

  const double computeMs = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - computeStart).count();
  gz::msgs::Double computeMsg;
  computeMsg.set_data(computeMs);
  this->computeTimePub.Publish(computeMsg);
}

GZ_ADD_PLUGIN(TendonArmController,
              gz::sim::System,
              gz::sim::ISystemConfigure,
              gz::sim::ISystemPreUpdate)
