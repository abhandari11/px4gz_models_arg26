#ifndef COSSERAT_ARM_CONTROLLER_COSSERAT_ARM_CONTROLLER_HH_
#define COSSERAT_ARM_CONTROLLER_COSSERAT_ARM_CONTROLLER_HH_

#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include <gz/sim/System.hh>
#include <gz/sim/Entity.hh>
#include <gz/sim/Model.hh>
#include <gz/transport/Node.hh>
#include <gz/msgs/vector3d.pb.h>
#include <gz/msgs/pose.pb.h>
#include <gz/msgs/double.pb.h>

#include "../common/CommandShaper.hpp"
#include "CosseratRodSolver.hpp"

namespace cosserat_arm_controller
{
/// Dynamic Cosserat-rod tendon controller for a 3-tendon soft continuum arm
/// (see models/continuum_arm_cosserat), replacing continuum_arm_pcc's
/// ad-hoc "split the total bend angle evenly across every joint" model with
/// a proper elasticity-based mechanics solve.
///
/// Pipeline, every physics tick:
///   0. Raw tendon commands (OnTendonCmd) do NOT drive the solver directly -
///      they only set the target of 3 CommandShaper filters (one per
///      tendon, see plugins/common/CommandShaper.hpp), which PreUpdate()
///      advances every tick. This gives the arm the same smooth,
///      zero-overshoot, ~2s-settling open-loop response as
///      continuum_arm_pcc, and turns the solver's input into a slowly-
///      ramping signal instead of an instant step.
///   1. At a throttled rate (solveRateHz, default ~50Hz - resolving on
///      every physics tick is unnecessary and was the source of "solver
///      chatter": re-solving to a hard instantaneous target every command
///      is what produced transient force spikes; a slowly-ramping input
///      converges each throttled solve in very few warm-started Newton
///      iterations instead), convert the CURRENT filtered tendon lengths
///      -> pulling tensions (a tendon can only pull, never push:
///      T_i = clamp(k_t * max(0, restLength - l_i), 0, maxTension)).
///   2. Solve the static Cosserat BVP (CosseratRodSolver::Solve) for the
///      resulting backbone shape (p(s), R(s)), warm-started from the
///      previous solve's converged base wrench.
///   3. Sample the solved backbone at the center of each of the N segments
///      (matching the SDF's own discretization) and convert each pair of
///      consecutive segment orientations into a 2-DOF (theta_x, theta_y)
///      relative-rotation target for that segment's universal joint.
///
/// Actuation itself reuses continuum_arm_pcc's proven mechanism (PD force
/// control directly against the real physics joints via the public
/// gz::sim::Joint API - see PreUpdate()) rather than kinematically
/// teleporting link poses: the segment links are still part of the SDF's
/// joint/kinematic tree, and DART's constraint solver only allows a joint's
/// child link to move within that joint's own DOF - an arbitrary
/// SetWorldPose() on an interior link would fight the constraint solver.
/// Driving the joint's own generalized coordinates via PD force is the
/// well-posed way to make the same continuous-mechanics result manifest as
/// real, physically simulated robot motion.
class CosseratArmController :
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

  /// gz-transport callback for the tendon length command topic. Only
  /// updates the CommandShapers' targets - see class header comment.
  private: void OnTendonCmd(const gz::msgs::Vector3d &_msg);

  /// Runs the Cosserat BVP solve for the given (already filtered) tendon
  /// lengths and updates jointTargets. Called from PreUpdate() at the
  /// throttled solveRateHz, NOT from OnTendonCmd.
  private: void RunSolve(double _l1, double _l2, double _l3);

  private: gz::sim::Model model{gz::sim::kNullEntity};
  private: std::vector<std::string> jointNames;
  private: std::vector<gz::sim::Entity> jointEntities;
  private: gz::sim::Entity tipLinkEntity{gz::sim::kNullEntity};

  private: std::unique_ptr<CosseratRodSolver> solver;
  private: double tendonStiffness{350.0};  // N/m of tendon shortening
  // N - safety clamp, empirically chosen (live-swept via gz_arm_harness,
  // see model.sdf's plugin comment) as the top of the LARGEST clean,
  // reversal-free branch of this rod's tension-vs-bend response starting
  // from tendonBaseline_=0. This rod's own Euler buckling load is <1N at
  // its E=5e5 Pa modulus, so tension-vs-bend is genuinely non-monotonic
  // (large-deflection elastica, not a solver bug) well below the servo's
  // own multi-hundred-N headroom, and a live sweep across E=5e5-4e6 Pa
  // confirmed a stiffer rod doesn't raise this ceiling (tension needed for
  // a given bend and the buckling threshold both scale with E together) -
  // max_tension is set by where the CURRENT branch's clean response tops
  // out, not by the servo, the rod's absolute strength, or its stiffness.
  private: double maxTension{5.25};
  // N, applied to ALL 3 tendons even when un-shortened (see OnTendonCmd()).
  // A live sweep found a nonzero baseline pretension SHRINKS the clean
  // (reversal-free) window vs. no baseline at all (0.18rad ceiling at
  // 5.5N baseline vs. 0.27rad ceiling at 0N baseline) - baseline pretension
  // on all 3 tendons simultaneously adds a constant net axial compressive
  // preload that eats into this rod's already-low (<1N) buckling margin
  // before any differential bending is even commanded. Kept at 0 for the
  // largest available clean branch; see maxTension's comment for where
  // that branch's own peak is.
  private: double tendonBaseline{0.0};
  private: double restLength{0.40};
  // If true, OnTendonCmd() calls CosseratRodSolver::SolveGuided() (with 3
  // discrete routing notches, chord-direction-consistent tendon loading)
  // instead of Solve()'s fixed-local-frame small-curvature approximation -
  // see CosseratRodSolver.hpp's SolveGuided() comment for why this raises
  // the clean (non-buckled) reachable bend range dramatically over the
  // small-curvature model this replaced.
  private: bool useGuidedTendons{false};

  // Retuned down from an earlier (6.0, 0.06) for the true hollow-tube
  // segment inertia, which is ~42% of a solid-rod segment's - the same
  // gains on the lower-inertia joints ran the local PD closed-loop near the
  // fixed-timestep physics integrator's stability boundary and visibly rang
  // instead of settling. Values below hold the same natural
  // frequency/damping ratio the old gains gave the old (solid-rod) inertia.
  private: double positionPGain{2.5};
  private: double positionDGain{0.05};
  // N*m - derived as maxTension * pitch_radius (16N * 0.015m = 0.24N*m):
  // the max torque the tendon-side actuation can realistically deliver to
  // any one joint, shared with continuum_arm_pcc's own effort_limit so
  // both models saturate the physics joints at the same servo-derived
  // torque ceiling (see continuum_arm_pcc/model.sdf's TendonArmController
  // plugin block).
  private: double effortLimit{0.24};

  private: std::string tendonTopic{"/continuum_arm_cosserat/tendon_cmd"};
  private: std::string tipPoseTopic{"/continuum_arm_cosserat/tip_pose"};

  /// Command-shaping filters (one per tendon) - see
  /// plugins/common/CommandShaper.hpp and class header comment. Shared
  /// defaults with continuum_arm_pcc's TendonArmController so both models'
  /// open-loop step responses are directly comparable.
  private: double commandFilterOmega{2.92};  // rad/s
  private: double maxTendonRateMPerS{0.06};  // m/s, 0 = unlimited
  private: continuum_arm_plugins::CommandShaper shaperL1;
  private: continuum_arm_plugins::CommandShaper shaperL2;
  private: continuum_arm_plugins::CommandShaper shaperL3;

  /// Throttled solve cadence - see class header comment item 1.
  private: double solveRateHz{50.0};
  private: double solveAccumulator{0.0};
  /// Re-solving a static (unchanged) input was found live to slowly drift
  /// the warm-started shooting solve off the correct straight-rod solution
  /// under gravity-only (zero-tension) loading - repeatedly re-seeding a
  /// Newton solve from its own immediately-prior (already-converged, but
  /// not bit-exact) iterate at 50Hz let a tiny per-solve residual asymmetry
  /// compound over hundreds of solves into tens of degrees of spurious
  /// rest-state droop; a single solve-once-per-command never re-fed its own
  /// output back in, so this never surfaced there. Only actually invoking
  /// RunSolve() when the filtered length has moved by more than this
  /// dead-band since the last solve removes the compounding path entirely
  /// while leaving an active ramp (which DOES keep moving every tick)
  /// solved at the full throttled rate.
  private: double solveDeadBandM{2e-4};
  private: bool haveSolved{false};
  private: double lastSolvedL1{0.0};
  private: double lastSolvedL2{0.0};
  private: double lastSolvedL3{0.0};

  private: gz::transport::Node node;
  private: gz::transport::Node::Publisher tipPosePub;
  /// Task-3 compute benchmark: wall-clock cost of RunSolve() (the BVP
  /// shooting solve), gz.msgs.Double milliseconds, published only when a
  /// solve actually runs (not every tick, unlike the other two models'
  /// trivial-cost PreUpdate() - a solve-only metric is the meaningful
  /// "cost of this model's mechanics" figure here).
  private: gz::transport::Node::Publisher computeTimePub;

  private: std::mutex mutex;
  /// Per-joint (theta_x, theta_y) targets, recomputed by RunSolve().
  private: std::vector<std::pair<double, double>> jointTargets;
};
}  // namespace cosserat_arm_controller

#endif
