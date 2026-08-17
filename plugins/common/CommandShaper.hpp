#ifndef CONTINUUM_ARM_PLUGINS_COMMON_COMMAND_SHAPER_HPP_
#define CONTINUUM_ARM_PLUGINS_COMMON_COMMAND_SHAPER_HPP_

#include <algorithm>

namespace continuum_arm_plugins
{

/// Header-only, per-scalar critically-damped 2nd-order command shaper,
/// shared by all three continuum-arm control plugins (tendon_arm_controller,
/// cosserat_arm_controller, rbsn_arm_controller) so raw tendon-length
/// commands never jump the plant's target instantaneously. One instance
/// per tendon (3 per model).
///
/// Integrates accel = omega^2*(target-value) - 2*omega*velocity at the
/// physics timestep (semi-implicit Euler). With damping ratio fixed at 1
/// (critically damped) the step response is mathematically monotonic -
/// zero overshoot by construction - with a smooth, bell-shaped velocity
/// profile (the "S-curve" shape requested for actuator command smoothing).
///
/// omega (rad/s) sets the response speed. For a critically-damped 2nd-order
/// system, x(t) = target - (target-x0)*(1+omega*t)*exp(-omega*t); solving
/// (1+omega*ts)*exp(-omega*ts) = 0.02 numerically for the 2% settling time
/// gives omega*ts ~= 5.83, i.e. omega ~= 5.83/ts. For the task's ts=2.0s
/// baseline this is omega ~= 2.92 rad/s - used as the default below, and
/// exposed as an SDF-tunable parameter per plugin so it can be verified/
/// adjusted live against the bench harness rather than trusted blindly.
///
/// An optional hard velocity clamp (maxRate, 0 = unlimited) sits on top,
/// modeling a physical actuator rate limit - set high enough that it stays
/// inactive at the calibration step (the filter alone should dominate
/// there) and only engages for larger commands.
class CommandShaper
{
  public: void Reset(double _value)
  {
    this->value_ = _value;
    this->velocity_ = 0.0;
    this->target_ = _value;
  }

  public: void SetTarget(double _target) { this->target_ = _target; }

  public: void SetOmega(double _omega) { this->omega_ = _omega; }

  public: void SetMaxRate(double _maxRate) { this->maxRate_ = _maxRate; }

  public: double Value() const { return this->value_; }

  public: double Velocity() const { return this->velocity_; }

  public: double Target() const { return this->target_; }

  /// Advances the filter by _dt seconds. Call once per physics tick.
  public: void Update(double _dt)
  {
    if (_dt <= 0.0)
      return;

    const double accel = this->omega_ * this->omega_ * (this->target_ - this->value_) -
                          2.0 * this->omega_ * this->velocity_;
    this->velocity_ += accel * _dt;
    if (this->maxRate_ > 0.0)
    {
      this->velocity_ =
          std::clamp(this->velocity_, -this->maxRate_, this->maxRate_);
    }
    this->value_ += this->velocity_ * _dt;
  }

  private: double value_{0.0};
  private: double velocity_{0.0};
  private: double target_{0.0};
  private: double omega_{2.92};   // rad/s, ~2% settle in ~2.0s (see above)
  private: double maxRate_{0.06}; // m/s, 0 = unlimited
};

}  // namespace continuum_arm_plugins

#endif
