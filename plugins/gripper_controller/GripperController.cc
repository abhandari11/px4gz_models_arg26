#include "GripperController.hh"

#include <sstream>

#include <gz/plugin/Register.hh>
#include <gz/common/Console.hh>

#include <gz/sim/Joint.hh>
#include <gz/sim/Link.hh>
#include <gz/sim/Util.hh>
#include <gz/sim/components/Name.hh>
#include <gz/sim/components/Static.hh>

#include <gz/msgs/stringmsg.pb.h>

#include <algorithm>

using namespace gripper_controller;

namespace
{
/// Finds a collision entity by name anywhere under a model, without needing
/// to know which link it belongs to.
gz::sim::Entity FindCollisionByName(
    gz::sim::Model &_model,
    gz::sim::EntityComponentManager &_ecm,
    const std::string &_name)
{
  for (const auto &link : _model.Links(_ecm))
  {
    gz::sim::Entity col = gz::sim::Link(link).CollisionByName(_ecm, _name);
    if (col != gz::sim::kNullEntity)
      return col;
  }
  return gz::sim::kNullEntity;
}
}  // namespace

//////////////////////////////////////////////////
void GripperController::Configure(
    const gz::sim::Entity &_entity,
    const std::shared_ptr<const sdf::Element> &_sdf,
    gz::sim::EntityComponentManager &_ecm,
    gz::sim::EventManager & /*_eventMgr*/)
{
  this->model = gz::sim::Model(_entity);
  if (!this->model.Valid(_ecm))
  {
    gzerr << "GripperController must be attached to a model entity."
          << std::endl;
    return;
  }

  if (_sdf->HasElement("vehicle_model_name"))
    this->vehicleModelName = _sdf->Get<std::string>("vehicle_model_name");

  std::string gripperBaseLinkName = "gripper_base";
  if (_sdf->HasElement("gripper_base_link"))
    gripperBaseLinkName = _sdf->Get<std::string>("gripper_base_link");
  this->gripperBaseLink = this->model.LinkByName(_ecm, gripperBaseLinkName);
  if (this->gripperBaseLink == gz::sim::kNullEntity)
  {
    gzerr << "GripperController: could not find gripper base link ["
          << gripperBaseLinkName << "]" << std::endl;
  }

  if (auto elem = _sdf->FindElement("finger_collision"))
  {
    while (elem)
    {
      this->fingerCollisionNames.push_back(elem->Get<std::string>());
      elem = elem->GetNextElement("finger_collision");
    }
  }
  for (const auto &name : this->fingerCollisionNames)
  {
    gz::sim::Entity col = FindCollisionByName(this->model, _ecm, name);
    if (col == gz::sim::kNullEntity)
    {
      gzerr << "GripperController: could not find finger collision ["
            << name << "]" << std::endl;
    }
    this->fingerCollisionEntities.push_back(col);
  }

  // Finger joints are actuated directly by this plugin (PD force control via
  // the public gz::sim::Joint API) rather than via the stock
  // JointPositionController system, which was found not to reliably move
  // these joints in every launch environment.
  if (auto elem = _sdf->FindElement("finger_joint"))
  {
    while (elem)
    {
      this->fingerJointNames.push_back(elem->Get<std::string>());
      elem = elem->GetNextElement("finger_joint");
    }
  }
  for (const auto &name : this->fingerJointNames)
  {
    gz::sim::Entity joint = this->model.JointByName(_ecm, name);
    if (joint == gz::sim::kNullEntity)
    {
      gzerr << "GripperController: could not find finger joint ["
            << name << "]" << std::endl;
    }
    else
    {
      gz::sim::Joint(joint).EnablePositionCheck(_ecm, true);
      gz::sim::Joint(joint).EnableVelocityCheck(_ecm, true);
    }
    this->fingerJointEntities.push_back(joint);
  }

  if (_sdf->HasElement("open_position"))
    this->openPosition = _sdf->Get<double>("open_position");
  if (_sdf->HasElement("closed_position"))
    this->closedPosition = _sdf->Get<double>("closed_position");
  if (_sdf->HasElement("grasp_radius"))
    this->graspRadius = _sdf->Get<double>("grasp_radius");
  if (_sdf->HasElement("grasp_center_offset"))
    this->graspCenterOffset = _sdf->Get<gz::math::Vector3d>("grasp_center_offset");
  if (_sdf->HasElement("finger_p_gain"))
    this->fingerPGain = _sdf->Get<double>("finger_p_gain");
  if (_sdf->HasElement("finger_d_gain"))
    this->fingerDGain = _sdf->Get<double>("finger_d_gain");
  if (_sdf->HasElement("finger_effort_limit"))
    this->fingerEffortLimit = _sdf->Get<double>("finger_effort_limit");
  if (_sdf->HasElement("command_topic"))
    this->commandTopic = _sdf->Get<std::string>("command_topic");
  if (_sdf->HasElement("state_topic"))
    this->stateTopic = _sdf->Get<std::string>("state_topic");

  if (auto elem = _sdf->FindElement("finger_contact_topic"))
  {
    while (elem)
    {
      this->fingerContactTopics.push_back(elem->Get<std::string>());
      elem = elem->GetNextElement("finger_contact_topic");
    }
  }
  this->lastContacts.resize(this->fingerContactTopics.size());

  this->statePub =
      this->node.Advertise<gz::msgs::StringMsg>(this->stateTopic);

  this->node.Subscribe(this->commandTopic, &GripperController::OnCommand,
                        this);

  for (std::size_t i = 0; i < this->fingerContactTopics.size(); ++i)
  {
    this->node.Subscribe(this->fingerContactTopics[i],
        std::function<void(const gz::msgs::Contacts &)>(
            [this, i](const gz::msgs::Contacts &_msg)
            {
              this->OnFingerContact(i, _msg);
            }));
  }

  gzmsg << "GripperController configured for model ["
        << this->vehicleModelName << "]: command topic ["
        << this->commandTopic << "], state topic [" << this->stateTopic
        << "]" << std::endl;
}

//////////////////////////////////////////////////
void GripperController::OnCommand(const gz::msgs::Boolean &_msg)
{
  std::lock_guard<std::mutex> lock(this->mutex);
  this->desiredClosed = _msg.data();
}

//////////////////////////////////////////////////
void GripperController::OnFingerContact(std::size_t _fingerIndex,
                                         const gz::msgs::Contacts &_msg)
{
  std::lock_guard<std::mutex> lock(this->mutex);
  if (_fingerIndex < this->lastContacts.size())
    this->lastContacts[_fingerIndex] = _msg;
}

//////////////////////////////////////////////////
bool GripperController::FingerTouches(
    gz::sim::EntityComponentManager &_ecm,
    std::size_t _fingerIndex,
    gz::sim::Entity _candidateModel)
{
  if (_fingerIndex >= this->fingerCollisionEntities.size())
    return false;

  gz::sim::Entity ownCollision = this->fingerCollisionEntities[_fingerIndex];

  gz::msgs::Contacts contacts;
  {
    std::lock_guard<std::mutex> lock(this->mutex);
    if (_fingerIndex >= this->lastContacts.size())
      return false;
    contacts = this->lastContacts[_fingerIndex];
  }

  for (int i = 0; i < contacts.contact_size(); ++i)
  {
    const auto &contact = contacts.contact(i);
    gz::sim::Entity c1 = contact.collision1().id();
    gz::sim::Entity c2 = contact.collision2().id();

    gz::sim::Entity other = gz::sim::kNullEntity;
    if (c1 == ownCollision)
      other = c2;
    else if (c2 == ownCollision)
      other = c1;
    else
      continue;

    if (gz::sim::topLevelModel(other, _ecm) == _candidateModel)
      return true;
  }
  return false;
}

//////////////////////////////////////////////////
void GripperController::Attach(
    gz::sim::EntityComponentManager &_ecm, gz::sim::Entity _object)
{
  auto gripperPose = gz::sim::Link(this->gripperBaseLink).WorldPose(_ecm);
  gz::sim::Model objectModel(_object);
  auto objectLink = objectModel.CanonicalLink(_ecm);
  auto objectPose = gz::sim::Link(objectLink).WorldPose(_ecm);
  if (!gripperPose.has_value() || !objectPose.has_value())
  {
    gzerr << "GripperController: could not resolve poses to attach object."
          << std::endl;
    return;
  }

  this->attachOffset = gripperPose->Inverse() * *objectPose;
  this->attached = true;
  this->attachedObject = _object;
  this->attachedObjectName = objectModel.Name(_ecm);

  // The pose-lock alone is what keeps the object "in the gripper" - we don't
  // need (or want) the fingers to also be in continuous physical contact
  // with it. Leaving collision on caused the fingers to keep reacting to the
  // object snapping back into place every tick (visible as random-looking
  // finger open/close jitter) and let the object's own dynamics fight the
  // pose override. Disabling collision while attached removes all of that;
  // Release() turns it back on so the object behaves normally once dropped.
  objectModel.SetCollisionEnabled(_ecm, false);

  gzmsg << "GripperController: grasped [" << this->attachedObjectName << "]"
        << std::endl;
}

//////////////////////////////////////////////////
void GripperController::Release(gz::sim::EntityComponentManager &_ecm)
{
  if (this->attached)
  {
    gzmsg << "GripperController: released [" << this->attachedObjectName
          << "]" << std::endl;
    gz::sim::Model(this->attachedObject).SetCollisionEnabled(_ecm, true);
  }
  this->attached = false;
  this->attachedObject = gz::sim::kNullEntity;
  this->attachedObjectName.clear();
}

//////////////////////////////////////////////////
void GripperController::DriveFingers(gz::sim::EntityComponentManager &_ecm)
{
  bool closed;
  {
    std::lock_guard<std::mutex> lock(this->mutex);
    closed = this->desiredClosed;
  }
  double target = closed ? this->closedPosition : this->openPosition;

  for (const auto &jointEntity : this->fingerJointEntities)
  {
    if (jointEntity == gz::sim::kNullEntity)
      continue;

    gz::sim::Joint joint(jointEntity);
    auto position = joint.Position(_ecm);
    auto velocity = joint.Velocity(_ecm);
    if (!position.has_value() || position->empty() ||
        !velocity.has_value() || velocity->empty())
    {
      continue;
    }

    double error = target - (*position)[0];
    double effort = this->fingerPGain * error - this->fingerDGain * (*velocity)[0];
    effort = std::clamp(effort, -this->fingerEffortLimit, this->fingerEffortLimit);
    joint.SetForce(_ecm, {effort});
  }
}

//////////////////////////////////////////////////
void GripperController::TryGrasp(gz::sim::EntityComponentManager &_ecm)
{
  auto gripperPose = gz::sim::Link(this->gripperBaseLink).WorldPose(_ecm);
  if (!gripperPose.has_value())
    return;

  gz::math::Vector3d graspCenter =
      gripperPose->Pos() + gripperPose->Rot() * this->graspCenterOffset;

  gz::sim::Entity worldEntity =
      gz::sim::worldEntity(this->model.Entity(), _ecm);
  gz::sim::World world(worldEntity);

  std::vector<gz::sim::Entity> validObjects;

  for (const auto &candidate : world.Models(_ecm))
  {
    gz::sim::Model candidateModel(candidate);

    // Ignore the UAV itself, and anything static (environment: ground,
    // walls, etc).
    if (candidateModel.Name(_ecm) == this->vehicleModelName)
      continue;
    if (candidateModel.Static(_ecm))
      continue;

    auto candidateLink = candidateModel.CanonicalLink(_ecm);
    auto candidatePose = gz::sim::Link(candidateLink).WorldPose(_ecm);
    if (!candidatePose.has_value())
      continue;

    double distance = (candidatePose->Pos() - graspCenter).Length();
    if (distance > this->graspRadius)
      continue;

    bool touchingAFinger = false;
    for (std::size_t i = 0; i < this->fingerCollisionEntities.size(); ++i)
    {
      if (this->FingerTouches(_ecm, i, candidate))
      {
        touchingAFinger = true;
        break;
      }
    }
    if (!touchingAFinger)
      continue;

    validObjects.push_back(candidate);
  }

  if (validObjects.size() == 1)
  {
    this->Attach(_ecm, validObjects.front());
  }
  else if (validObjects.size() > 1)
  {
    gzwarn << "GripperController: grasp ambiguous (" << validObjects.size()
           << " candidates), not attaching." << std::endl;
  }
  // validObjects.empty(): nothing in reach yet - not an error, just keep
  // waiting (PreUpdate retries this every tick while closed and unattached,
  // since the fingers take time to physically swing shut after the command
  // arrives - checking only once, right when the command is received, would
  // almost always miss objects that aren't already touching at that instant).
}

//////////////////////////////////////////////////
void GripperController::PreUpdate(
    const gz::sim::UpdateInfo & /*_info*/, gz::sim::EntityComponentManager &_ecm)
{
  bool closed;
  {
    std::lock_guard<std::mutex> lock(this->mutex);
    closed = this->desiredClosed;
  }

  this->DriveFingers(_ecm);

  int closedInt = closed ? 1 : 0;
  if (closedInt != this->lastCommandedClosed)
  {
    if (!closed)
      this->Release(_ecm);
    this->lastCommandedClosed = closedInt;
  }

  // Keep retrying every tick while closed and not yet holding anything: the
  // fingers need time to physically swing shut after the command arrives, so
  // checking only once (right on the open->closed edge) would almost always
  // miss objects that weren't already touching at that exact instant.
  if (closed && !this->attached)
    this->TryGrasp(_ecm);

  if (this->attached)
  {
    auto gripperPose = gz::sim::Link(this->gripperBaseLink).WorldPose(_ecm);
    if (gripperPose.has_value())
    {
      gz::sim::Model attachedModel(this->attachedObject);
      attachedModel.SetWorldPoseCmd(_ecm, *gripperPose * this->attachOffset);
      // Collision is off while attached, but the object is still a normal
      // dynamic body underneath - zero its velocity each tick too so it
      // doesn't quietly build up momentum (e.g. from gravity) that would
      // suddenly show up as an unexpected jump/throw on release.
      gz::sim::Entity objectLink = attachedModel.CanonicalLink(_ecm);
      gz::sim::Link(objectLink).SetLinearVelocity(_ecm, {0, 0, 0});
      gz::sim::Link(objectLink).SetAngularVelocity(_ecm, {0, 0, 0});
    }
  }
}

//////////////////////////////////////////////////
void GripperController::PostUpdate(
    const gz::sim::UpdateInfo & /*_info*/,
    const gz::sim::EntityComponentManager & /*_ecm*/)
{
  bool closed;
  {
    std::lock_guard<std::mutex> lock(this->mutex);
    closed = this->desiredClosed;
  }

  std::ostringstream json;
  json << "{\"closed\": " << (closed ? "true" : "false")
       << ", \"attached\": " << (this->attached ? "true" : "false")
       << ", \"object\": \"" << this->attachedObjectName << "\"}";

  gz::msgs::StringMsg stateMsg;
  stateMsg.set_data(json.str());
  this->statePub.Publish(stateMsg);
}

GZ_ADD_PLUGIN(GripperController,
              gz::sim::System,
              gz::sim::ISystemConfigure,
              gz::sim::ISystemPreUpdate,
              gz::sim::ISystemPostUpdate)
