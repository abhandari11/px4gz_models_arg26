#ifndef TENDON_ARM_CONTROLLER_TENDON_ARM_CONTROLLER_HH_
#define TENDON_ARM_CONTROLLER_TENDON_ARM_CONTROLLER_HH_

#include <mutex>
#include <string>
#include <vector>

#include <gz/sim/System.hh>
#include <gz/sim/Entity.hh>
#include <gz/sim/Model.hh>
#include <gz/transport/Node.hh>
#include <gz/msgs/vector3d.pb.h>
#include <gz/msgs/double.pb.h>

#include "../common/CommandShaper.hpp"

namespace tendon_arm_controller
{
/// Constant-curvature tendon-space joint-mimic controller for a 3-tendon
/// soft continuum arm (see models/continuum_arm_pcc).
///
/// Subscribes to a commanded 3-tendon length vector (l1,l2,l3). Raw
/// commands are NOT applied instantly - each of the 3 tendon channels is
/// first shaped by a critically-damped CommandShaper (see
/// plugins/common/CommandShaper.hpp), giving the arm the smooth,
/// zero-overshoot, ~2s-settling actuator response the open-loop baseline
/// targets. The shaped (filtered) lengths are what actually feed the
/// existing tendon-to-curvature Jacobian: maps to global bend angles
/// (theta_x, theta_y), distributes theta_x/N and theta_y/N evenly across
/// every joint listed in <joint_name>, and drives each joint's two axes
/// toward that target with a PD force loop. The per-joint passive
/// spring/damper declared in the SDF <axis><dynamics> block still acts
/// underneath this - the PD loop is the "actuation" layered on top of the
/// silicone's own passive restoring force, not a replacement for it.
class TendonArmController :
  public gz::sim::System,
  public gz::sim::ISystemConfigure,
  public gz::sim::ISystemPreUpdate
{
  public: void Configure(
              const gz::sim::Entity &_entity,
              const std::shared_ptr<const sdf::Element> &_sdf,
              gz::sim::EntityComponentManager &_ecm,
              gz::sim::EventManager &_eventMgr) override;

  public: void PreUpdate(
              const gz::sim::UpdateInfo &_info,
              gz::sim::EntityComponentManager &_ecm) override;

  /// gz-transport callback for the tendon length command topic.
  private: void OnTendonCmd(const gz::msgs::Vector3d &_msg);

  /// Recomputes thetaXTarget/thetaYTarget from the latest tendon command.
  /// Called with the mutex held.
  private: void UpdateTargets();

  private: gz::sim::Model model{gz::sim::kNullEntity};
  private: std::vector<std::string> jointNames;
  private: std::vector<gz::sim::Entity> jointEntities;

  private: double pitchRadius{0.015};
  private: double restLength{0.40};
  private: double positionPGain{0.15};
  private: double positionDGain{0.004};
  private: double effortLimit{0.3};
  private: std::string tendonTopic{"/continuum_arm_pcc/tendon_cmd"};
  private: std::string anglesTopic{"/continuum_arm_pcc/bend_angles"};

  /// Command-shaping filters (one per tendon) - see
  /// plugins/common/CommandShaper.hpp. Raw commands set their target only;
  /// PreUpdate() advances them every tick and feeds the FILTERED length
  /// into UpdateTargets(), not the raw command.
  private: double commandFilterOmega{2.92};  // rad/s
  private: double maxTendonRateMPerS{0.06};  // m/s, 0 = unlimited
  private: continuum_arm_plugins::CommandShaper shaperL1;
  private: continuum_arm_plugins::CommandShaper shaperL2;
  private: continuum_arm_plugins::CommandShaper shaperL3;
  private: bool shapersInitialized{false};

  private: gz::transport::Node node;
  private: gz::transport::Node::Publisher anglesPub;
  /// Task-3 compute benchmark: wall-clock cost of this plugin's own
  /// PreUpdate() work (gz.msgs.Double, milliseconds), published every
  /// tick so the eval harness can aggregate mean/peak per model. PCC's
  /// own PD/Jacobian evaluation is trivial (no numerical solve), so this
  /// is mainly a baseline for comparison against Cosserat's much costlier
  /// BVP solve.
  private: gz::transport::Node::Publisher computeTimePub;

  private: std::mutex mutex;
  /// Latest FILTERED tendon lengths (l1, l2, l3) - written by PreUpdate()
  /// each tick from the CommandShapers, then consumed by UpdateTargets().
  /// Defaults to restLength each (i.e. zero deflection) until the first
  /// command arrives.
  private: double l1{0.0};
  private: double l2{0.0};
  private: double l3{0.0};
  /// Raw (unfiltered) commanded tendon lengths - the CommandShapers' target.
  private: double l1Cmd{0.0};
  private: double l2Cmd{0.0};
  private: double l3Cmd{0.0};
  private: bool haveCommand{false};

  /// Global bend angles (rad), and their per-joint share, recomputed
  /// whenever a new tendon command arrives.
  private: double thetaX{0.0};
  private: double thetaY{0.0};
  private: double thetaXPerJoint{0.0};
  private: double thetaYPerJoint{0.0};
};
}  // namespace tendon_arm_controller

#endif
