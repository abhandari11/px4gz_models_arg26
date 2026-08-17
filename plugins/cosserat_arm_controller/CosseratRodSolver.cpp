#include "CosseratRodSolver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace cosserat_arm_controller
{

//////////////////////////////////////////////////
double Vec3::Norm() const { return std::sqrt(x * x + y * y + z * z); }

//////////////////////////////////////////////////
Mat3 Mat3::operator*(const Mat3 &o) const
{
  Mat3 r;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      r.m[i][j] = m[i][0] * o.m[0][j] + m[i][1] * o.m[1][j] + m[i][2] * o.m[2][j];
  return r;
}

Mat3 Mat3::Transpose() const
{
  Mat3 r;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      r.m[i][j] = m[j][i];
  return r;
}

Mat3 Mat3::operator+(const Mat3 &o) const
{
  Mat3 r;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      r.m[i][j] = m[i][j] + o.m[i][j];
  return r;
}

Mat3 Mat3::operator*(double s) const
{
  Mat3 r;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      r.m[i][j] = m[i][j] * s;
  return r;
}

void Mat3::Orthonormalize()
{
  // Modified Gram-Schmidt on the 3 column vectors.
  Vec3 c0(m[0][0], m[1][0], m[2][0]);
  Vec3 c1(m[0][1], m[1][1], m[2][1]);
  Vec3 c2(m[0][2], m[1][2], m[2][2]);

  double n0 = c0.Norm();
  if (n0 < 1e-12) n0 = 1.0;
  c0 = c0 * (1.0 / n0);

  c1 = c1 - c0 * c0.Dot(c1);
  double n1 = c1.Norm();
  if (n1 < 1e-12) n1 = 1.0;
  c1 = c1 * (1.0 / n1);

  c2 = c0.Cross(c1);  // enforce right-handedness exactly

  m[0][0] = c0.x; m[0][1] = c1.x; m[0][2] = c2.x;
  m[1][0] = c0.y; m[1][1] = c1.y; m[1][2] = c2.y;
  m[2][0] = c0.z; m[2][1] = c1.z; m[2][2] = c2.z;
}

//////////////////////////////////////////////////
Mat3 Hat(const Vec3 &u)
{
  Mat3 r;
  r.m[0][0] = 0;    r.m[0][1] = -u.z; r.m[0][2] = u.y;
  r.m[1][0] = u.z;  r.m[1][1] = 0;    r.m[1][2] = -u.x;
  r.m[2][0] = -u.y; r.m[2][1] = u.x;  r.m[2][2] = 0;
  return r;
}

//////////////////////////////////////////////////
CosseratRodSolver::CosseratRodSolver(const Params &_params) : params_(_params)
{
  // Hollow tube (annulus) envelope cross-section: A = pi*(ro^2-ri^2),
  // I = (pi/4)*(ro^4-ri^4) about any diametral axis, J = 2I. solidityFillFactor
  // then scales these down for the constitutive (Kse_/Kbt_) computation only,
  // modeling the wall as a sparse strut lattice rather than solid material -
  // see Params::solidityFillFactor.
  const double ro = params_.outerRadius;
  const double ri = params_.innerRadius;
  const double A = M_PI * (ro * ro - ri * ri);
  const double I = (M_PI / 4.0) * (ro * ro * ro * ro - ri * ri * ri * ri);
  const double J = 2.0 * I;                     // polar moment

  const double eta = params_.solidityFillFactor;
  const double Aeff = eta * A;
  const double Ieff = eta * I;
  const double Jeff = eta * J;

  this->Kse_ = Vec3(params_.shearModulus * Aeff, params_.shearModulus * Aeff,
                     params_.youngsModulus * Aeff);
  this->Kbt_ = Vec3(params_.youngsModulus * Ieff, params_.youngsModulus * Ieff,
                     params_.shearModulus * Jeff);

  // Mass-per-length intentionally uses the FULL envelope area A, not Aeff:
  // `density` is already an effective/homogenized value over the full
  // envelope (see Params::density's own comment), so it must not be scaled
  // by the fill factor a second time.
  const double massPerLength = params_.density * A;
  this->gravityForcePerLength_ = Vec3(0, 0, -9.81) * massPerLength;
}

//////////////////////////////////////////////////
CosseratRodSolver::State CosseratRodSolver::Derivative(const State &_y) const
{
  const Vec3 RtN = _y.R.Transpose() * _y.n;
  const Vec3 v(this->vStar_.x + RtN.x / this->Kse_.x,
               this->vStar_.y + RtN.y / this->Kse_.y,
               this->vStar_.z + RtN.z / this->Kse_.z);

  const Vec3 RtM = _y.R.Transpose() * _y.m;
  const Vec3 u(this->uStar_.x + RtM.x / this->Kbt_.x,
               this->uStar_.y + RtM.y / this->Kbt_.y,
               this->uStar_.z + RtM.z / this->Kbt_.z);

  State d;
  d.p = _y.R * v;
  d.R = _y.R * Hat(u);

  // Gravity is always a true distributed (per-unit-length) load.
  // tendonForcePerLength_/tendonMomentPerLength_ are an ADDITIONAL
  // distributed load, populated by Solve() only when
  // Params::distributeTendonLoadContinuously is set (zero otherwise, in
  // which case tendon loading stays a concentrated boundary condition
  // instead - see class header comment).
  d.n = -this->gravityForcePerLength_ - _y.R * this->tendonForcePerLength_;
  d.m = -(d.p.Cross(_y.n)) - _y.R * this->tendonMomentPerLength_;
  return d;
}

//////////////////////////////////////////////////
CosseratRodSolver::State CosseratRodSolver::Integrate(
    const Vec3 &_n0, const Vec3 &_m0,
    const std::vector<NotchLoad> &_notches)
{
  State y;
  y.p = Vec3(0, 0, 0);
  y.R = Mat3::Identity();
  y.n = _n0;
  y.m = _m0;

  const int steps = params_.integrationSteps;
  const double h = params_.length / steps;

  this->trajectory_.clear();
  this->sValues_.clear();
  this->trajectory_.reserve(steps + 1);
  this->sValues_.reserve(steps + 1);
  this->trajectory_.push_back(y);
  this->sValues_.push_back(0.0);

  std::size_t nextNotch = 0;

  for (int i = 0; i < steps; ++i)
  {
    const State k1 = Derivative(y);

    State y2;
    y2.p = y.p + k1.p * (h / 2.0);
    y2.R = y.R + k1.R * (h / 2.0);
    y2.n = y.n + k1.n * (h / 2.0);
    y2.m = y.m + k1.m * (h / 2.0);
    const State k2 = Derivative(y2);

    State y3;
    y3.p = y.p + k2.p * (h / 2.0);
    y3.R = y.R + k2.R * (h / 2.0);
    y3.n = y.n + k2.n * (h / 2.0);
    y3.m = y.m + k2.m * (h / 2.0);
    const State k3 = Derivative(y3);

    State y4;
    y4.p = y.p + k3.p * h;
    y4.R = y.R + k3.R * h;
    y4.n = y.n + k3.n * h;
    y4.m = y.m + k3.m * h;
    const State k4 = Derivative(y4);

    y.p = y.p + (k1.p + k2.p * 2.0 + k3.p * 2.0 + k4.p) * (h / 6.0);
    y.R = y.R + (k1.R + k2.R * 2.0 + k3.R * 2.0 + k4.R) * (h / 6.0);
    y.n = y.n + (k1.n + k2.n * 2.0 + k3.n * 2.0 + k4.n) * (h / 6.0);
    y.m = y.m + (k1.m + k2.m * 2.0 + k3.m * 2.0 + k4.m) * (h / 6.0);
    y.R.Orthonormalize();

    const double sNow = (i + 1) * h;

    // Interior tendon-guide jumps (see class header comment / Solve()):
    // each notch this step has just reached or passed gets its own
    // force/moment subtracted from the running internal wrench, exactly
    // the same natural-BC form the tip itself uses (n(L) = R(L)*forceLocal),
    // just applied partway along the domain instead of only at s=L.
    while (nextNotch < _notches.size() && sNow >= _notches[nextNotch].s - 1e-9)
    {
      y.n = y.n - y.R * _notches[nextNotch].forceLocal;
      y.m = y.m - y.R * _notches[nextNotch].momentLocal;
      ++nextNotch;
    }

    this->trajectory_.push_back(y);
    this->sValues_.push_back(sNow);
  }

  return y;
}

//////////////////////////////////////////////////
bool CosseratRodSolver::Solve(const std::vector<double> &_tendonTensions,
                               const Vec3 &_tipForce, const Vec3 &_tipMoment)
{
  // Combined (summed over tendons) per-notch-point LOCAL force/moment,
  // under the small-curvature "tendon tangent ~= local tangent (+z)"
  // approximation (see class header comment). Constant per Solve() call
  // (depends only on tension/geometry); tendonForceLocal/tendonMomentLocal
  // is each tendon's FULL pull, then split evenly across
  // params_.tendonNotchFractions below (see class header comment for why
  // -- neither a single concentrated tip load nor a fully continuous
  // smeared load, but a small number of discrete guide points matching
  // the physical hardware's routing notches).
  Vec3 tendonForceLocal(0, 0, 0);
  Vec3 tendonMomentLocal(0, 0, 0);
  const int numTendons = static_cast<int>(_tendonTensions.size());
  const double d = params_.tendonPitchRadius;
  for (int i = 0; i < numTendons; ++i)
  {
    // Negative sign: this solver's local frame is related to the SDF
    // model's actual body frame by a fixed 180-degree flip about the
    // shared local X axis (see CosseratArmController.cc's OnTendonCmd()
    // comment on the matching thetaY correction) - a solver-frame CCW
    // tendon layout maps to a CW one in the real model frame. Using -theta
    // here keeps commanded tendon index i's bend landing at the SAME
    // labeled real-world direction (2*pi*i/numTendons) the test/caller
    // expects, without which tendons 1 and 2 (not on the shared X axis)
    // would swap which real-world direction they bend toward.
    const double theta = -2.0 * M_PI * i / numTendons;
    const Vec3 ri(d * std::cos(theta), d * std::sin(theta), 0.0);
    const double T = _tendonTensions[i];
    // Tendon pulls toward the base along its own (~local +z) direction at
    // each guide/anchor point: force contribution -T*z. Moment
    // contribution from the offset is +T*(ri x z), not the
    // naively-expected -T*(ri x z): the shooting method's warm-started
    // root for this rod (well past its own <1N Euler buckling load at any
    // working tension) lands on the elastica branch with this
    // empirically-verified sign - see the tension-sweep scratchpad
    // tooling referenced in model.sdf, and CosseratArmController.cc's
    // OnTendonCmd() comment on the paired thetaY correction this must
    // stay consistent with. (Verified against this multi-notch loading
    // too, not just the original single-tip-BC version.)
    tendonForceLocal = tendonForceLocal - Vec3(0, 0, T);
    tendonMomentLocal = tendonMomentLocal + ri.Cross(Vec3(0, 0, T));
  }

  Vec3 notchForceLocal(0, 0, 0);
  Vec3 notchMomentLocal(0, 0, 0);
  std::vector<NotchLoad> interiorNotches;
  if (params_.distributeTendonLoadContinuously)
  {
    // Fold the WHOLE combined tendon force/moment into the ODE as a
    // uniform per-unit-length term (see Derivative()) instead of any
    // point load - the tip boundary condition below then has nothing
    // tendon-related left to balance (notchForceLocal/notchMomentLocal
    // stay zero, interiorNotches stays empty).
    // NOTE: tendonMomentLocal's sign (+T*(ri x z), see the per-tendon
    // loop's own comment) was empirically tuned for the CONCENTRATED tip
    // boundary condition specifically -- test showed distributed mode
    // reverses the bend direction from what that sign gives, so flip it
    // here rather than at the source (keeps the concentrated-mode path,
    // which is verified correct, untouched).
    this->tendonForcePerLength_ = tendonForceLocal * (1.0 / params_.length);
    this->tendonMomentPerLength_ = tendonMomentLocal * (-1.0 / params_.length);
  }
  else
  {
    this->tendonForcePerLength_ = Vec3(0, 0, 0);
    this->tendonMomentPerLength_ = Vec3(0, 0, 0);

    // Split evenly across the notch points; the LAST point is the tip,
    // still handled via the shooting method's own natural BC, while every
    // earlier point is an explicit jump applied during Integrate().
    const int numPoints = std::max<int>(1, static_cast<int>(params_.tendonNotchFractions.size()));
    notchForceLocal = tendonForceLocal * (1.0 / numPoints);
    notchMomentLocal = tendonMomentLocal * (1.0 / numPoints);
    for (int k = 0; k + 1 < numPoints; ++k)
      interiorNotches.push_back(
          {params_.tendonNotchFractions[k] * params_.length, notchForceLocal, notchMomentLocal});
  }

  // Unknowns x = [n0(3), m0(3)]. Natural BC at s=L: the rod's own internal
  // wrench there must balance the externally applied tip wrench (this
  // notch's own share of the tendon anchor load -- zero when the tendon
  // load is folded into the ODE distributedly instead, see above --
  // rotated into world frame via the trial's own tip orientation, plus
  // any additional _tipForce/_tipMoment) -
  // E(x) = [n(L) - R(L)*notchForceLocal - tipForce;
  //         m(L) - R(L)*notchMomentLocal - tipMoment].
  double x[6] = {this->n0_.x, this->n0_.y, this->n0_.z,
                 this->m0_.x, this->m0_.y, this->m0_.z};

  auto residual = [&](const double _x[6]) -> std::array<double, 6>
  {
    const Vec3 n0(_x[0], _x[1], _x[2]);
    const Vec3 m0(_x[3], _x[4], _x[5]);
    const State tip = Integrate(n0, m0, interiorNotches);
    const Vec3 worldTendonForce = tip.R * notchForceLocal;
    const Vec3 worldTendonMoment = tip.R * notchMomentLocal;
    const Vec3 fErr = tip.n - worldTendonForce - _tipForce;
    const Vec3 mErr = tip.m - worldTendonMoment - _tipMoment;
    return {fErr.x, fErr.y, fErr.z, mErr.x, mErr.y, mErr.z};
  };

  bool converged = false;
  std::array<double, 6> bestResidual{};
  double bestNorm = -1.0;
  double bestX[6] = {0, 0, 0, 0, 0, 0};

  for (int iter = 0; iter < params_.maxNewtonIterations; ++iter)
  {
    const std::array<double, 6> E = residual(x);
    double normE = 0.0;
    for (double e : E) normE += e * e;
    normE = std::sqrt(normE);

    if (bestNorm < 0.0 || normE < bestNorm)
    {
      bestNorm = normE;
      bestResidual = E;
      for (int k = 0; k < 6; ++k) bestX[k] = x[k];
    }

    if (normE < params_.newtonTolerance)
    {
      converged = true;
      break;
    }

    // Finite-difference 6x6 Jacobian dE/dx.
    double J[6][6];
    for (int k = 0; k < 6; ++k)
    {
      double xPerturbed[6];
      for (int j = 0; j < 6; ++j) xPerturbed[j] = x[j];
      const double eps = (k < 3) ? params_.fdEpsilonForce : params_.fdEpsilonMoment;
      xPerturbed[k] += eps;
      const std::array<double, 6> Ep = residual(xPerturbed);
      for (int row = 0; row < 6; ++row)
        J[row][k] = (Ep[row] - E[row]) / eps;
    }

    // Solve J * dx = -E via Gaussian elimination with partial pivoting.
    double A[6][7];
    for (int row = 0; row < 6; ++row)
    {
      for (int col = 0; col < 6; ++col) A[row][col] = J[row][col];
      A[row][6] = -E[row];
    }
    bool singular = false;
    for (int col = 0; col < 6 && !singular; ++col)
    {
      int pivot = col;
      double pivotMag = std::fabs(A[col][col]);
      for (int row = col + 1; row < 6; ++row)
      {
        if (std::fabs(A[row][col]) > pivotMag)
        {
          pivot = row;
          pivotMag = std::fabs(A[row][col]);
        }
      }
      if (pivotMag < 1e-14)
      {
        singular = true;
        break;
      }
      if (pivot != col)
        for (int k = 0; k < 7; ++k) std::swap(A[col][k], A[pivot][k]);

      for (int row = 0; row < 6; ++row)
      {
        if (row == col) continue;
        const double factor = A[row][col] / A[col][col];
        for (int k = col; k < 7; ++k) A[row][k] -= factor * A[col][k];
      }
    }

    if (singular)
      break;  // Fall back to best iterate found so far.

    double dx[6];
    for (int row = 0; row < 6; ++row) dx[row] = A[row][6] / A[row][row];

    // Backtracking line search: a full (or fixed-damping) Newton step can
    // overshoot badly when starting far from the root (e.g. cold start at
    // x=0 with a large commanded tension, observed empirically to diverge
    // with a plain fixed-damping step). Halving the step until the
    // residual norm actually decreases is the standard, simple fix.
    double stepScale = 1.0;
    double xTrial[6];
    double normTrial = normE;
    for (int backtrack = 0; backtrack < 12; ++backtrack)
    {
      for (int k = 0; k < 6; ++k) xTrial[k] = x[k] + stepScale * dx[k];
      const std::array<double, 6> ETrial = residual(xTrial);
      double n2 = 0.0;
      for (double e : ETrial) n2 += e * e;
      normTrial = std::sqrt(n2);
      if (normTrial < normE || stepScale < 1e-4)
        break;
      stepScale *= 0.5;
    }
    for (int k = 0; k < 6; ++k) x[k] = xTrial[k];
  }

  // Re-integrate once more with the best x found so that this->trajectory_
  // (used by SampleSegmentCenters()) reflects the returned solution even
  // if the loop above exited on the residual-evaluation call rather than
  // ending with a fresh Integrate().
  Integrate(Vec3(bestX[0], bestX[1], bestX[2]),
            Vec3(bestX[3], bestX[4], bestX[5]),
            interiorNotches);

  this->n0_ = Vec3(bestX[0], bestX[1], bestX[2]);
  this->m0_ = Vec3(bestX[3], bestX[4], bestX[5]);

  return converged;
}

//////////////////////////////////////////////////
std::vector<CosseratRodSolver::Sample> CosseratRodSolver::SampleSegmentCenters(
    int _n) const
{
  std::vector<Sample> out;
  if (this->trajectory_.empty() || _n <= 0)
    return out;
  out.reserve(_n);

  const double segLen = params_.length / _n;
  std::size_t searchStart = 0;
  for (int i = 0; i < _n; ++i)
  {
    const double sTarget = (i + 0.5) * segLen;

    // Find the bracketing pair in the dense RK4 trajectory and linearly
    // interpolate position/orientation (orientation interpolation is a
    // simple per-column lerp + re-orthonormalize - the dense trajectory
    // has integrationSteps/N sub-steps per segment, typically >= 5, so the
    // angular change between adjacent samples is small enough that this is
    // a good approximation without needing full SLERP).
    std::size_t j = searchStart;
    while (j + 1 < this->sValues_.size() && this->sValues_[j + 1] < sTarget)
      ++j;
    searchStart = j;

    const std::size_t j1 = std::min(j + 1, this->sValues_.size() - 1);
    const double s0 = this->sValues_[j];
    const double s1 = this->sValues_[j1];
    const double t = (s1 > s0) ? (sTarget - s0) / (s1 - s0) : 0.0;

    Sample sample;
    sample.s = sTarget;
    sample.p = this->trajectory_[j].p * (1.0 - t) + this->trajectory_[j1].p * t;
    sample.R = this->trajectory_[j].R * (1.0 - t) + this->trajectory_[j1].R * t;
    sample.R.Orthonormalize();
    out.push_back(sample);
  }
  return out;
}

//////////////////////////////////////////////////
std::vector<CosseratRodSolver::Sample> CosseratRodSolver::NotchBackboneSamples() const
{
  std::vector<Sample> out;
  if (this->trajectory_.empty())
    return out;

  std::size_t searchStart = 0;
  for (double frac : params_.tendonNotchFractions)
  {
    const double sTarget = frac * params_.length;
    std::size_t j = searchStart;
    while (j + 1 < this->sValues_.size() && this->sValues_[j + 1] < sTarget)
      ++j;
    searchStart = j;

    const std::size_t j1 = std::min(j + 1, this->sValues_.size() - 1);
    const double s0 = this->sValues_[j];
    const double s1 = this->sValues_[j1];
    const double t = (s1 > s0) ? (sTarget - s0) / (s1 - s0) : 0.0;

    Sample sample;
    sample.s = sTarget;
    sample.p = this->trajectory_[j].p * (1.0 - t) + this->trajectory_[j1].p * t;
    sample.R = this->trajectory_[j].R * (1.0 - t) + this->trajectory_[j1].R * t;
    sample.R.Orthonormalize();
    out.push_back(sample);
  }
  return out;
}

//////////////////////////////////////////////////
bool CosseratRodSolver::SolveGuided(const std::vector<double> &_tendonTensions,
                                     int _maxOuterIterations, double _outerTolerance)
{
  const int numTendons = static_cast<int>(_tendonTensions.size());
  const double d = params_.tendonPitchRadius;
  const int numNotches = static_cast<int>(params_.tendonNotchFractions.size());

  std::vector<double> thetas(numTendons);
  std::vector<Vec3> ri(numTendons);
  for (int i = 0; i < numTendons; ++i)
  {
    // Same routing-angle convention as Solve() (see its per-tendon loop
    // comment on the -theta sign).
    thetas[i] = -2.0 * M_PI * i / numTendons;
    ri[i] = Vec3(d * std::cos(thetas[i]), d * std::sin(thetas[i]), 0.0);
  }

  // Warm-start from the PREVIOUS call's converged per-notch forces (see
  // this method's header comment) when available and dimensionally
  // compatible; otherwise fall back to Solve()'s small-curvature
  // approximation (uniform local +z pull) as the cold-start seed. Reused
  // as-is even though the tension usually changed since the previous call
  // - the outer loop's own iteration corrects both direction and
  // magnitude, and starting from an already-curved direction converges in
  // a handful of passes instead of the ~40-70+ a from-scratch axial guess
  // needs.
  std::vector<std::vector<NotchLoad>> perTendonNotches;
  const bool haveWarmStart =
      static_cast<int>(this->guidedNotchesWarmStart_.size()) == numTendons &&
      numTendons > 0 &&
      static_cast<int>(this->guidedNotchesWarmStart_[0].size()) == numNotches;
  if (haveWarmStart)
  {
    perTendonNotches = this->guidedNotchesWarmStart_;
  }
  else
  {
    perTendonNotches.resize(numTendons);
    for (int i = 0; i < numTendons; ++i)
    {
      const double T = _tendonTensions[i];
      const Vec3 f(0, 0, -T);
      const Vec3 mo = ri[i].Cross(Vec3(0, 0, T));
      perTendonNotches[i].resize(numNotches);
      for (int k = 0; k < numNotches; ++k)
        perTendonNotches[i][k] = {params_.tendonNotchFractions[k] * params_.length, f, mo};
    }
  }

  bool innerConverged = false;
  bool outerConverged = false;
  for (int outer = 0; outer < _maxOuterIterations; ++outer)
  {
    std::vector<NotchLoad> interior;
    Vec3 tipForceLocal(0, 0, 0), tipMomentLocal(0, 0, 0);
    for (int k = 0; k < numNotches; ++k)
    {
      Vec3 f(0, 0, 0), mo(0, 0, 0);
      for (int i = 0; i < numTendons; ++i)
      {
        f = f + perTendonNotches[i][k].forceLocal;
        mo = mo + perTendonNotches[i][k].momentLocal;
      }
      if (k + 1 < numNotches)
        interior.push_back({params_.tendonNotchFractions[k] * params_.length, f, mo});
      else
      {
        tipForceLocal = f;
        tipMomentLocal = mo;
      }
    }

    // Inner Newton BVP solve - same machinery as Solve(), just against the
    // (possibly notch-specific) `interior` list and explicit tip
    // force/moment instead of one shared value.
    double x[6] = {this->n0_.x, this->n0_.y, this->n0_.z,
                   this->m0_.x, this->m0_.y, this->m0_.z};
    auto residual = [&](const double _x[6]) -> std::array<double, 6>
    {
      const Vec3 n0(_x[0], _x[1], _x[2]);
      const Vec3 m0(_x[3], _x[4], _x[5]);
      const State tip = Integrate(n0, m0, interior);
      const Vec3 fErr = tip.n - tip.R * tipForceLocal;
      const Vec3 mErr = tip.m - tip.R * tipMomentLocal;
      return {fErr.x, fErr.y, fErr.z, mErr.x, mErr.y, mErr.z};
    };

    bool converged = false;
    double bestNorm = -1.0;
    double bestX[6] = {0, 0, 0, 0, 0, 0};

    for (int iter = 0; iter < params_.maxNewtonIterations; ++iter)
    {
      const std::array<double, 6> E = residual(x);
      double normE = 0.0;
      for (double e : E) normE += e * e;
      normE = std::sqrt(normE);

      if (bestNorm < 0.0 || normE < bestNorm)
      {
        bestNorm = normE;
        for (int k = 0; k < 6; ++k) bestX[k] = x[k];
      }
      if (normE < params_.newtonTolerance)
      {
        converged = true;
        break;
      }

      double J[6][6];
      for (int k = 0; k < 6; ++k)
      {
        double xPerturbed[6];
        for (int j = 0; j < 6; ++j) xPerturbed[j] = x[j];
        const double eps = (k < 3) ? params_.fdEpsilonForce : params_.fdEpsilonMoment;
        xPerturbed[k] += eps;
        const std::array<double, 6> Ep = residual(xPerturbed);
        for (int row = 0; row < 6; ++row)
          J[row][k] = (Ep[row] - E[row]) / eps;
      }

      double A[6][7];
      for (int row = 0; row < 6; ++row)
      {
        for (int col = 0; col < 6; ++col) A[row][col] = J[row][col];
        A[row][6] = -E[row];
      }
      bool singular = false;
      for (int col = 0; col < 6 && !singular; ++col)
      {
        int pivot = col;
        double pivotMag = std::fabs(A[col][col]);
        for (int row = col + 1; row < 6; ++row)
        {
          if (std::fabs(A[row][col]) > pivotMag)
          {
            pivot = row;
            pivotMag = std::fabs(A[row][col]);
          }
        }
        if (pivotMag < 1e-14)
        {
          singular = true;
          break;
        }
        if (pivot != col)
          for (int k = 0; k < 7; ++k) std::swap(A[col][k], A[pivot][k]);
        for (int row = 0; row < 6; ++row)
        {
          if (row == col) continue;
          const double factor = A[row][col] / A[col][col];
          for (int k = col; k < 7; ++k) A[row][k] -= factor * A[col][k];
        }
      }
      if (singular)
        break;

      double dx[6];
      for (int row = 0; row < 6; ++row) dx[row] = A[row][6] / A[row][row];

      double stepScale = 1.0;
      double xTrial[6];
      double normTrial = normE;
      for (int backtrack = 0; backtrack < 12; ++backtrack)
      {
        for (int k = 0; k < 6; ++k) xTrial[k] = x[k] + stepScale * dx[k];
        const std::array<double, 6> ETrial = residual(xTrial);
        double n2 = 0.0;
        for (double e : ETrial) n2 += e * e;
        normTrial = std::sqrt(n2);
        if (normTrial < normE || stepScale < 1e-4)
          break;
        stepScale *= 0.5;
      }
      for (int k = 0; k < 6; ++k) x[k] = xTrial[k];
    }

    innerConverged = converged;
    this->n0_ = Vec3(bestX[0], bestX[1], bestX[2]);
    this->m0_ = Vec3(bestX[3], bestX[4], bestX[5]);
    Integrate(this->n0_, this->m0_, interior);

    // Recompute each tendon's per-notch loads from the ACTUAL solved shape:
    // at each notch, a frictionless routing point feels T toward whichever
    // neighboring anchor (previous notch, or the base entry point for the
    // first notch) plus T toward the next neighbor (or nothing, for the
    // LAST notch, where the tendon terminates instead of continuing on).
    // This is what makes notches mutually reinforcing instead of each
    // independently assuming local +z: a notch downstream of a real bend
    // now pulls along the ACTUAL chord to its neighbors, not a fixed
    // straight-ahead guess - see class header / SolveGuided()'s own
    // declaration comment for why this fixes the buckling collapse the
    // fixed-local-frame approximation produced.
    std::vector<Sample> backbone = NotchBackboneSamples();
    double maxChange = 0.0;
    for (int i = 0; i < numTendons; ++i)
    {
      std::vector<Vec3> anchors(numNotches);
      for (int k = 0; k < numNotches; ++k)
        anchors[k] = backbone[k].p + backbone[k].R * ri[i];

      const double T = _tendonTensions[i];
      Vec3 prevAnchor = ri[i];  // base entry point (R(0)=I, p(0)=0)
      for (int k = 0; k < numNotches; ++k)
      {
        Vec3 dirIn = prevAnchor - anchors[k];
        const double nIn = dirIn.Norm();
        if (nIn > 1e-9) dirIn = dirIn * (1.0 / nIn);

        Vec3 dirOut(0, 0, 0);
        if (k + 1 < numNotches)
        {
          dirOut = anchors[k + 1] - anchors[k];
          const double nOut = dirOut.Norm();
          if (nOut > 1e-9) dirOut = dirOut * (1.0 / nOut);
        }

        const Vec3 worldForce = (dirIn + dirOut) * T;
        const Vec3 localForce = backbone[k].R.Transpose() * worldForce;
        const Vec3 localMoment = ri[i].Cross(localForce);

        const double change = (localForce - perTendonNotches[i][k].forceLocal).Norm();
        maxChange = std::max(maxChange, change);

        // Damped update (outer fixed-point relaxation) - a full-strength
        // update was found to oscillate between two chord estimates on
        // successive passes for highly-bent configurations instead of
        // settling; 0.25 converges reliably within a few hundred passes.
        const double relax = 0.25;
        perTendonNotches[i][k].forceLocal =
            perTendonNotches[i][k].forceLocal * (1.0 - relax) + localForce * relax;
        perTendonNotches[i][k].momentLocal =
            perTendonNotches[i][k].momentLocal * (1.0 - relax) + localMoment * relax;

        prevAnchor = anchors[k];
      }
    }

    if (maxChange < _outerTolerance)
    {
      outerConverged = true;
      break;
    }
  }

  // Store for the next call's warm start regardless of whether this call
  // fully converged within budget - even a partially-relaxed estimate is a
  // better starting point than the cold small-curvature guess once the
  // rod has any real curvature.
  this->guidedNotchesWarmStart_ = perTendonNotches;

  return outerConverged && innerConverged;
}

}  // namespace cosserat_arm_controller
