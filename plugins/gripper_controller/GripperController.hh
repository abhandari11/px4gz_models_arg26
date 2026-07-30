#ifndef GRIPPER_CONTROLLER_GRIPPER_CONTROLLER_HH_
#define GRIPPER_CONTROLLER_GRIPPER_CONTROLLER_HH_

#include <mutex>
#include <string>
#include <vector>

#include <gz/sim/System.hh>
#include <gz/sim/Entity.hh>
#include <gz/sim/Model.hh>
#include <gz/sim/World.hh>
#include <gz/math/Pose3.hh>
#include <gz/transport/Node.hh>
#include <gz/msgs/boolean.pb.h>
#include <gz/msgs/contacts.pb.h>

namespace gripper_controller
{
/// Translates a single open/close command into finger motion and dynamic
/// grasp/release logic for a procedural gripper mounted on a model.
///
/// All link/joint/topic names are read from SDF plugin parameters (not
/// hardcoded) so a future CAD gripper can be swapped in by editing the SDF
/// block only, without touching this class.
///
/// The "grasp" is implemented as a kinematic pose-lock (the grasped object's
/// world pose is written directly, tracking the gripper base each tick)
/// rather than a true dynamically-created physics joint: gz-sim8's stock
/// DetachableJoint system only attaches a statically-declared child at load
/// time, and there is no verified public API for creating a joint between
/// two independently-spawned models at runtime. The pose-lock gives the same
/// externally-visible behavior (rigid attachment, clean release) using only
/// documented, public ECM APIs.
class GripperController :
  public gz::sim::System,
  public gz::sim::ISystemConfigure,
  public gz::sim::ISystemPreUpdate,
  public gz::sim::ISystemPostUpdate
{
  public: void Configure(
              const gz::sim::Entity &_entity,
              const std::shared_ptr<const sdf::Element> &_sdf,
              gz::sim::EntityComponentManager &_ecm,
              gz::sim::EventManager &_eventMgr) override;

  public: void PreUpdate(
              const gz::sim::UpdateInfo &_info,
              gz::sim::EntityComponentManager &_ecm) override;

  public: void PostUpdate(
              const gz::sim::UpdateInfo &_info,
              const gz::sim::EntityComponentManager &_ecm) override;

  /// Drives all finger joints toward the currently-commanded target angle
  /// with a simple PD force controller. Implemented directly against the
  /// public gz::sim::Joint API (SetForce/Position/Velocity) rather than the
  /// stock JointPositionController system, since that stock plugin was found
  /// not to actuate these joints reliably in every launch environment.
  private: void DriveFingers(gz::sim::EntityComponentManager &_ecm);

  /// Search the world for a single graspable object inside the grasp volume
  /// and in contact with a finger, and attach it if exactly one is found.
  private: void TryGrasp(gz::sim::EntityComponentManager &_ecm);

  /// Attach (pose-lock) the given object entity to the gripper.
  private: void Attach(gz::sim::EntityComponentManager &_ecm,
                        gz::sim::Entity _object);

  /// Release the currently-attached object, if any.
  private: void Release(gz::sim::EntityComponentManager &_ecm);

  /// gz-transport callback for the open/close command topic.
  private: void OnCommand(const gz::msgs::Boolean &_msg);

  /// gz-transport callback for a finger's contact-sensor topic.
  private: void OnFingerContact(std::size_t _fingerIndex,
                                 const gz::msgs::Contacts &_msg);

  /// True if any recent contact on the given finger involves the given
  /// collision entity's owning top-level model.
  private: bool FingerTouches(gz::sim::EntityComponentManager &_ecm,
                               std::size_t _fingerIndex,
                               gz::sim::Entity _candidateModel);

  private: gz::sim::Model model{gz::sim::kNullEntity};
  private: gz::sim::Entity gripperBaseLink{gz::sim::kNullEntity};
  private: std::vector<std::string> fingerJointNames;
  private: std::vector<gz::sim::Entity> fingerJointEntities;
  private: std::vector<std::string> fingerCollisionNames;
  private: std::vector<gz::sim::Entity> fingerCollisionEntities;

  private: std::string vehicleModelName{"hillstar"};
  private: double openPosition{0.0};
  private: double closedPosition{0.9};
  private: double graspRadius{0.08};
  /// Offset (in gripper_base_link's own local frame) from that link's origin
  /// to the actual "palm center" where the fingers do their work - needed
  /// because gripper_base_link is usually the vehicle's main body link, whose
  /// own origin can be a long way from the fingertips. Without this, the
  /// grasp-radius check measures from the body center instead of where
  /// contact actually happens, and never finds anything in range.
  private: gz::math::Vector3d graspCenterOffset{0, 0, 0};
  private: double fingerPGain{1.5};
  private: double fingerDGain{0.02};
  private: double fingerEffortLimit{0.5};
  private: std::string commandTopic{"/hillstar/gripper_cmd"};
  private: std::string stateTopic{"/hillstar/gripper_state"};
  private: std::vector<std::string> fingerContactTopics;

  private: gz::transport::Node node;
  private: gz::transport::Node::Publisher statePub;

  private: std::mutex mutex;
  private: bool desiredClosed{false};
  /// Tri-state edge tracking: -1 = not yet initialized.
  private: int lastCommandedClosed{-1};
  private: std::vector<gz::msgs::Contacts> lastContacts;

  private: bool attached{false};
  private: gz::sim::Entity attachedObject{gz::sim::kNullEntity};
  private: std::string attachedObjectName;
  private: gz::math::Pose3d attachOffset;
};
}  // namespace gripper_controller

#endif
