#ifndef COSSERAT_ARM_CONTROLLER_COSSERAT_ROD_SOLVER_HPP_
#define COSSERAT_ARM_CONTROLLER_COSSERAT_ROD_SOLVER_HPP_

#include <array>
#include <vector>

namespace cosserat_arm_controller
{

/// Minimal, dependency-free 3-vector/3x3-matrix pair used by the solver.
/// Kept independent of gz::math so this file is a self-contained numerical
/// component (statics of a tendon-driven Cosserat rod) that could be
/// dropped into any host, tested standalone, etc. - the gz-sim plugin layer
/// (CosseratArmController) is the only place that talks to Gazebo types.
struct Vec3
{
  double x{0.0}, y{0.0}, z{0.0};

  Vec3() = default;
  Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

  Vec3 operator+(const Vec3 &o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3 &o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator-() const { return {-x, -y, -z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

  double Dot(const Vec3 &o) const { return x * o.x + y * o.y + z * o.z; }
  Vec3 Cross(const Vec3 &o) const
  {
    return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
  }
  double Norm() const;
};

/// Row-major 3x3 matrix (used for the backbone orientation R(s)).
struct Mat3
{
  double m[3][3]{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

  static Mat3 Identity() { return Mat3(); }

  Vec3 operator*(const Vec3 &v) const
  {
    return {
        m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
        m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
        m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
    };
  }
  Mat3 operator*(const Mat3 &o) const;
  Mat3 Transpose() const;
  Mat3 operator+(const Mat3 &o) const;
  Mat3 operator*(double s) const;

  /// Gram-Schmidt re-orthonormalization - RK4-integrated rotation matrices
  /// drift off SO(3) over many steps; call this after every integration
  /// step to keep R a valid rotation (standard practice for this kind of
  /// matrix-ODE rod integration).
  void Orthonormalize();
};

/// Skew-symmetric (hat) map: so(3) vector -> 3x3 matrix, s.t. Hat(u)*v == u x v.
Mat3 Hat(const Vec3 &u);

/// Static, tendon-driven Cosserat rod BVP solver.
///
/// Solves the classic Cosserat rod statics ODEs (Rucker & Webster 2011,
/// global/world-frame n,m convention):
///   p'(s) = R(s) v(s)
///   R'(s) = R(s) Hat(u(s))
///   n'(s) = -f_ext(s)
///   m'(s) = -(p'(s) x n(s)) - l_ext(s)
/// with the linear constitutive law
///   v(s) = v* + Kse^-1 R(s)^T n(s)
///   u(s) = u* + Kbt^-1 R(s)^T m(s)
/// via a shooting method: guess the base wrench (n(0), m(0)), integrate
/// with RK4 to s=L, and Newton-Raphson-correct the guess against the tip
/// wrench boundary condition (finite-difference Jacobian, since an
/// analytical sensitivity equation isn't worth the complexity for a 6x6
/// root-find run once per tendon command).
///
/// Tendon modeling: each tendon threads through tendonNotchFractions.size()
/// discrete guide/anchor points along the rod (default: thirds, matching
/// the physical hardware's 3 internal routing notches), terminating at the
/// LAST point (s=L, the tip -- still the shooting method's own natural BC).
/// The tendon's total tension is split evenly across these points (T/N
/// each), applied at every point via the same small-curvature local
/// force/moment formula (see Solve()), rather than either of the two
/// extremes already tried: (a) the FULL tension concentrated only at the
/// tip, which put this soft (E=5e5 Pa) hollow tube's ENTIRE 0.4m length
/// under the tendon's full compressive load from the very first mm of
/// pull -- multiple times past its own <1N Euler buckling threshold at
/// any operationally useful bend, producing a converged-but-genuinely-
/// buckled multi-inflection-point backbone shape instead of a clean
/// single curve (verified live: per-joint curvature magnitude oscillates
/// rather than tapering monotonically, reproducibly, at every tension
/// past ~0.7N); or (b) a continuous body-frame load smeared uniformly
/// along the whole rod, which an earlier version of this solver used and
/// found even worse (peaked around 10-20N and then fell). Splitting
/// across a SMALL number of discrete points, applied during integration
/// as a state jump in n(s)/m(s) at each point's own s-location (see
/// Integrate()), keeps the total base reaction equal to the tendon's
/// actual tension (no artificial force multiplication) while building the
/// compressive load up gradually toward the base instead of applying it
/// at full magnitude over the rod's entire unsupported length -- directly
/// matching how the real tendon is mechanically supported at each notch,
/// not just anchored once at the tip.
///
/// This solver uses the small-curvature approximation that each tendon's
/// local tangent at each guide/anchor cross-section stays close to the
/// backbone's own local tangent there (valid since d=0.015m is small
/// compared to the rod's bend radius over the joint-limited curvatures
/// this arm operates at) - i.e. tendon i's local direction at each point
/// is taken as approximately +z rather than solving the fully implicit
/// direction-dependent formulation. That reduces each point's own share
/// of the tendon's contribution to a constant local-frame force
/// (-(Ti/N)*z) and moment (-(Ti/N)*(ri x z)), applied at the tip via the
/// natural BC n(L)=R(L)*forceLocal, m(L)=R(L)*momentLocal (see Solve()),
/// and at every other point via an explicit jump mid-integration.
class CosseratRodSolver
{
  public: struct Params
  {
    double length{0.40};             // m, total arc length L
    // Rod cross-section is a hollow tube (outer/inner radius), matching the
    // real arm: Ø40mm OD / Ø30mm ID lattice-TPU tube with 3 tendons routed
    // through the inner wall. area/Ixx/Iyy/J below are the true annulus
    // values, not a solid-rod approximation - see constructor.
    double outerRadius{0.020};       // m, 40mm OD
    double innerRadius{0.015};       // m, 30mm ID
    // Grid-pattern mesh fill factor (0-1): the real tube's wall is a sparse
    // strut lattice, not a continuous solid annulus wall - photographed
    // hardware shows a rigid-plastic crossed-strut lattice with roughly
    // 80% open area. Applied to A/Ixx/Iyy/J (below, in the constructor) but
    // deliberately NOT to mass-per-length: `density` above is already an
    // effective/homogenized value over the FULL envelope area (see its own
    // comment), so scaling mass by this factor too would double-count the
    // lattice's openness. Default 1.0 keeps the old solid-annulus behavior.
    double solidityFillFactor{1.0};
    double youngsModulus{5.0e5};     // Pa - see solidityFillFactor above:
                                      // this is the LATTICE STRUT MATERIAL's
                                      // own modulus (not a pre-softened
                                      // "effective" value) - the openness is
                                      // captured geometrically via
                                      // solidityFillFactor instead. A rigid
                                      // printed plastic (PLA/PETG-class,
                                      // E_solid ~ 1.5-2.5 GPa) at a sparse
                                      // ~15-20% fill lands in the 1e7-1e8 Pa
                                      // range for the resulting EI - see
                                      // continuum_arm_cosserat/model.sdf's
                                      // plugin block for the exact derived
                                      // value (fit to match the PCC model's
                                      // own K_theta~=1.611 N*m/rad baseline).
    double shearModulus{1.6e5};      // Pa
    double density{780.0};           // kg/m^3, effective bulk density of the
                                      // lattice-infill TPU (~65% of solid
                                      // TPU's ~1200 kg/m^3), sized so the
                                      // full arm+base assembly totals 250g
                                      // - see gen_cosserat_arm.py
    double tendonPitchRadius{0.015}; // m, 3 tendons on the 30mm-diameter
                                      // inner wall, 120deg apart
    // Discrete tendon guide/anchor points along the rod, as fractions of
    // L, matching the physical hardware's routing notches. Must be
    // strictly increasing and end at 1.0 (the tip, still the shooting
    // method's own boundary condition) - see class header comment.
    //
    // Defaults to a single point (the tip only, i.e. mathematically
    // identical to this solver's original tip-only-BC behavior). Solve()
    // -- the small-curvature approximation, every notch using a fixed
    // local-+z tendon direction regardless of actual rod shape -- made
    // curvature-profile waviness WORSE, not better, when this was
    // live-tested at {1/3, 2/3, 1}: each notch's force/moment was computed
    // independently in that same fixed local frame, so as soon as the rod
    // had ANY curvature the 3 points' contributions stopped reinforcing
    // one consistent pull and started geometrically interfering with each
    // other. SolveGuided() is the properly self-consistent fix for that
    // (tracks each notch's REAL chord direction to its neighbors, from the
    // rod's own current shape, instead of assuming local +z) - set this to
    // {1/3, 2/3, 1} (or similar) when calling SolveGuided(); leave it at
    // the {1.0} default for plain Solve() calls, which still only supports
    // the fixed-local-frame approximation.
    std::vector<double> tendonNotchFractions{1.0};
    // If true, ignore tendonNotchFractions and instead apply the tendons'
    // combined force/moment as a UNIFORM per-unit-length term folded
    // directly into the n'(s)/m'(s) ODE (see Derivative()), with the tip
    // boundary condition reduced to just any extra _tipForce/_tipMoment
    // (zero by default) - i.e. no single point (tip or otherwise) carries
    // a concentrated load; every cross-section only ever reacts the
    // PORTION of the total tendon force applied at s' > s (a classical
    // "column under distributed axial load" rather than "column under a
    // tip load", which has a higher critical buckling load for the same
    // total force since the peak compressive stress region is effectively
    // shorter). Simpler to implement than the discrete notch machinery
    // above (one ODE term, no jump/breakpoint bookkeeping) - see class
    // header comment for why this is being tried as an alternative to the
    // notch-point approach instead of built on top of it.
    bool distributeTendonLoadContinuously{false};
    int integrationSteps{150};       // RK4 steps over [0, L]
    int maxNewtonIterations{40};
    double newtonTolerance{1e-4};    // on ||tip wrench residual||
    double fdEpsilonForce{1e-4};     // finite-difference step, N
    double fdEpsilonMoment{1e-5};    // finite-difference step, N*m
  };

  public: struct Sample
  {
    double s{0.0};
    Vec3 p;
    Mat3 R;
  };

  public: explicit CosseratRodSolver(const Params &_params);

  /// Solves the BVP for the given tendon tensions (N, one per tendon,
  /// pulling only - caller is expected to clamp negative/slack tensions to
  /// 0 before calling), applied as a concentrated tip force/moment (see
  /// class header comment), plus any additional external tip wrench.
  /// Tendons are assumed evenly spaced (2*pi/numTendons apart) starting at
  /// angle 0. Returns true if Newton-Raphson converged within tolerance;
  /// on non-convergence the best (lowest-residual) iterate found is still
  /// stored/sampleable, so callers can degrade gracefully rather than
  /// freezing.
  public: bool Solve(const std::vector<double> &_tendonTensions,
                      const Vec3 &_tipForce = Vec3(0, 0, 0),
                      const Vec3 &_tipMoment = Vec3(0, 0, 0));

  /// Samples the last-solved backbone at `_n` points evenly spaced over
  /// (0, L] - i.e. at the centers of `_n` equal-length segments, matching
  /// how the N discretized links in the SDF are laid out.
  public: std::vector<Sample> SampleSegmentCenters(int _n) const;

  /// Last converged (or best-effort) base wrench - exposed so callers can
  /// warm-start the next Solve() call with it (successive tendon commands
  /// are usually close to each other, so warm-starting cuts Newton
  /// iterations dramatically).
  public: Vec3 BaseForce() const { return this->n0_; }
  public: Vec3 BaseMoment() const { return this->m0_; }

  /// Self-consistent guided-tendon solve: each tendon's per-notch force is
  /// recomputed every outer iteration from the ACTUAL chord direction
  /// between that notch's own routing-offset position and its neighbors'
  /// (previous notch or base anchor; next notch, or nothing for the last
  /// notch, where the tendon terminates), instead of Solve()'s fixed
  /// local-+z small-curvature approximation. This is the "properly
  /// self-consistent guided-tendon model" the Params::tendonNotchFractions
  /// comment flagged as not-yet-attempted: the earlier discrete-notch
  /// attempt used the SAME local-frame approximation at every notch
  /// regardless of how much the rod had already bent by that point, which
  /// is exactly what made multiple notches interfere with each other
  /// instead of reinforcing one consistent pull (see that comment, and
  /// CosseratArmController.cc's plugin config comment, for the buckling
  /// investigation this fixes: tracking each notch's real chord direction
  /// makes the tendon's force at each guide progressively MORE lateral and
  /// LESS axially-compressive as the rod bends, a self-reinforcing
  /// mechanism that live/offline-solver testing showed avoids the
  /// small-curvature model's multi-lobe buckling collapse for far higher
  /// bend angles than the un-guided model can reach cleanly). An outer
  /// fixed-point loop (damped, relax=0.25) alternates: (1) inner Newton
  /// BVP solve with the current per-notch forces, (2) recompute those
  /// forces from the resulting shape's actual notch positions - until the
  /// per-notch force stops changing. Returns true only if BOTH the outer
  /// loop and the final inner Newton solve converged.
  ///
  /// Warm-starts the outer loop itself from the PREVIOUS call's converged
  /// per-notch forces (member state, like n0_/m0_ already do for the inner
  /// Newton solve) rather than re-deriving the small-curvature seed every
  /// call - successive real commands are usually close to each other, and
  /// without this a cold re-seed needs tens of outer passes (~40ms) EVERY
  /// call, which is too slow for a realtime control loop publishing
  /// commands faster than that (found via live closed-loop testing: the
  /// solver falls behind and the plant looks like it never reaches its
  /// target, even though the same tension converges to the correct shape
  /// in an offline single-shot test). Falls back to the small-curvature
  /// seed only on the very first call (or after a tendon-count change).
  public: bool SolveGuided(const std::vector<double> &_tendonTensions,
                            int _maxOuterIterations = 300,
                            double _outerTolerance = 1e-3);

  private: struct State
  {
    Vec3 p;
    Mat3 R;
    Vec3 n;
    Vec3 m;
  };

  /// One interior notch's own local force/moment jump - possibly distinct
  /// per notch (SolveGuided()) or the same shared value at every notch
  /// (Solve()'s legacy small-curvature approximation).
  private: struct NotchLoad
  {
    double s{0.0};
    Vec3 forceLocal;
    Vec3 momentLocal;
  };

  /// Right-hand side of the ODE system. Gravity is always a true
  /// distributed load; tendonForcePerLength_/tendonMomentPerLength_ are
  /// an ADDITIONAL distributed load, nonzero only when
  /// Params::distributeTendonLoadContinuously is set (see Solve(), which
  /// populates them each call) - otherwise tendon loading stays a
  /// concentrated boundary condition and these are left at zero.
  private: State Derivative(const State &_y) const;

  /// Integrates from s=0 (given base wrench) to s=L with fixed-step RK4,
  /// storing the full trajectory in this->trajectory_ for later sampling.
  /// Applies a state jump (n -= R*forceLocal, m -= R*momentLocal) at every
  /// notch in _notches crossed during integration - the interior tendon
  /// guide points (excludes the tip, which is the shooting method's own
  /// boundary condition instead - see Solve()). Each notch may carry its
  /// own distinct force/moment (SolveGuided() gives every notch a
  /// different, chord-direction-derived value; Solve() gives them all the
  /// same legacy small-curvature-approximation value).
  private: State Integrate(const Vec3 &_n0, const Vec3 &_m0,
                            const std::vector<NotchLoad> &_notches);

  /// Samples the backbone (p(s), R(s)) from the last-integrated trajectory
  /// at every s in tendonNotchFractions*L (see SolveGuided(), which offsets
  /// these by each tendon's own routing radius to get that tendon's actual
  /// world-frame anchor path).
  private: std::vector<Sample> NotchBackboneSamples() const;

  private: Params params_;

  // Precomputed (constructor-time) constitutive/loading constants.
  private: Vec3 Kse_;                   // diag(GA, GA, EA)
  private: Vec3 Kbt_;                   // diag(E*Ixx, E*Iyy, G*J)
  private: Vec3 gravityForcePerLength_; // rho*A*g, world frame
  private: const Vec3 vStar_{0, 0, 1};  // stress-free linear strain (unit tangent)
  private: const Vec3 uStar_{0, 0, 0};  // stress-free angular strain (straight)

  // Set fresh by Solve() every call (from that call's tendon tensions);
  // read by Derivative(). Local (cross-section) frame, like
  // tendonForceLocal/tendonMomentLocal in Solve() - stays zero unless
  // Params::distributeTendonLoadContinuously is set.
  private: Vec3 tendonForcePerLength_{0, 0, 0};
  private: Vec3 tendonMomentPerLength_{0, 0, 0};

  private: Vec3 n0_{0, 0, 0};
  private: Vec3 m0_{0, 0, 0};
  private: std::vector<State> trajectory_;
  private: std::vector<double> sValues_;

  // Warm-start state for SolveGuided() - see its declaration comment.
  private: std::vector<std::vector<NotchLoad>> guidedNotchesWarmStart_;
};

}  // namespace cosserat_arm_controller

#endif
