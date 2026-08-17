#ifndef RBSN_ARM_CONTROLLER_RBSN_ARM_CONTROLLER_HH_
#define RBSN_ARM_CONTROLLER_RBSN_ARM_CONTROLLER_HH_

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

namespace rbsn_arm_controller
{
/// Non-linear Rigid-Body Spring Network (RBSN) continuum-arm controller
/// (see models/continuum_arm_rbsn) - the third of three architecturally
/// distinct mechanics models in this repo (kinematic PD-tracking for PCC,
/// quasi-static BVP+PD-tracking for Cosserat). Unlike those two, RBSN runs
/// GENUINE forward dynamics: no solved/kinematic target is ever position-
/// tracked. Instead, every physics tick, tendon tension is converted
/// directly into a generalized actuation TORQUE via the tendon-routing
/// Jacobian's transpose, and each of the 15 rigid segment joints carries
/// its own non-linear elastic restoring torque
/// tau_spring(theta) = k1*theta + k2*theta^3 plus linear viscous damping
/// tau_damp = c_theta*theta_dot - both applied directly as
/// joint.SetForce(tau_actuation - tau_spring - tau_damp), letting the
/// gz/DART physics engine integrate the resulting motion. There is no
/// position-error PD term anywhere in this plugin.
///
/// Pipeline, every physics tick:
///   1. Raw tendon length commands (OnTendonCmd) only set 3 CommandShaper
///      targets (see plugins/common/CommandShaper.hpp) - same command-
///      shaping filter continuum_arm_pcc/continuum_arm_cosserat use, for a
///      directly-comparable open-loop baseline across all three models.
///   2. The filtered (slowly-ramping) tendon lengths convert to pulling
///      tensions via the same tendon-stiffness spring-clutch model
///      continuum_arm_cosserat uses: T_i = clamp(tendon_baseline +
///      tendon_stiffness*max(0, rest_length - l_i), 0, max_tension).
///   3. Tensions convert to a total generalized actuation torque
///      (tau_x, tau_y) via the tendon-routing Jacobian's transpose (see
///      RbsnArmController.cc's header comment for the derivation), split
///      evenly across the 15 joints (matching PCC/Cosserat's own even-
///      split actuation-distribution convention).
///   4. Each joint's own non-linear spring + damper reacts against that
///      actuation torque; net force applied via joint.SetForce(). The
///      resulting equilibrium (actuation torque balances tau_spring at
///      some non-zero bend) emerges from real dynamics, not a solved
///      target - this is what makes RBSN a genuinely different model from
///      the other two, not a reskin of either.
class RbsnArmController :
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

  /// gz-transport callback for the tendon length command topic. Only
  /// updates the CommandShapers' targets - see class header comment.
  private: void OnTendonCmd(const gz::msgs::Vector3d &_msg);

  private: gz::sim::Model model{gz::sim::kNullEntity};
  private: std::vector<std::string> jointNames;
  private: std::vector<gz::sim::Entity> jointEntities;

  private: double pitchRadius{0.015};
  private: double restLength{0.40};

  // Tendon tension model - shares continuum_arm_cosserat's already-
  // validated actuation assumptions (same physical tendon/material class);
  // only k1/k2/cTheta below are RBSN-specific free parameters per the
  // task's own scope.
  private: double tendonStiffness{3333.0}; // N/m of tendon shortening
  private: double maxTension{100.0};       // N
  private: double tendonBaseline{0.0};     // N, pretension on all 3 tendons

  // Non-linear per-joint spring: tau_spring(theta) = k1*theta + k2*theta^3.
  private: double k1{1.6};    // N*m/rad
  private: double k2{0.0};    // N*m/rad^3
  private: double cTheta{0.05}; // N*m*s/rad, viscous joint damping

  private: double effortLimit{3.0}; // N*m, safety clamp on net per-joint torque

  private: std::string tendonTopic{"/continuum_arm_rbsn/tendon_cmd"};

  // Command-shaping filters (one per tendon) - see
  // plugins/common/CommandShaper.hpp. Shared defaults with the other two
  // models' identical filter so all three models' open-loop step
  // responses are directly comparable.
  private: double commandFilterOmega{2.92}; // rad/s
  private: double maxTendonRateMPerS{0.06}; // m/s, 0 = unlimited
  private: continuum_arm_plugins::CommandShaper shaperL1;
  private: continuum_arm_plugins::CommandShaper shaperL2;
  private: continuum_arm_plugins::CommandShaper shaperL3;

  private: gz::transport::Node node;
  /// Task-3 compute benchmark: wall-clock cost of this plugin's own
  /// PreUpdate() work (gz.msgs.Double, milliseconds), published every
  /// tick - RBSN's per-tick cost is a direct force evaluation (no
  /// numerical solve at all), so this is expected to land near PCC's own
  /// trivial cost and well below Cosserat's BVP solve cost.
  private: gz::transport::Node::Publisher computeTimePub;

  private: std::mutex mutex;
};
}  // namespace rbsn_arm_controller

#endif
