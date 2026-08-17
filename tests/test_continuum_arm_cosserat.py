"""
Functional/verification test for the continuum_arm_cosserat model +
CosseratArmController plugin (dynamic Cosserat rod solver replacing
continuum_arm_pcc's piecewise-constant-curvature model).

Launches `gz sim` headless with worlds/continuum_arm_cosserat_test.sdf (the
arm rigidly mounted to the world, clear of the ground), publishes tendon
commands on /continuum_arm_cosserat/tendon_cmd, and:

  1. Direction check (3 cases): shortening each tendon alone should deflect
     the tip TOWARD that tendon's own radial direction - same as a
     bimetallic strip curling toward its shorter fiber: the shortened side
     becomes the concave (inside) side of the bend, so the tip ends up on
     that same side. Read from the plugin's own /continuum_arm_cosserat/
     tip_pose (published every tick from the actual, physically simulated
     tip link pose - not the solver's raw/internal output), so this is a
     genuine "did gazebo physics actually move the way the mechanics call
     for" check, not a re-test of the solver's math in isolation.

  2. Continuity check: ramps a single tendon's commanded shortening up in
     several steps and logs the tip position at each step, verifying the
     deflection magnitude increases monotonically (a discontinuous jump or
     a reversal would indicate a BVP shooting failure/wrong-branch solution
     for some step - see CosseratRodSolver.hpp's header comment on the
     small-curvature tendon approximation's validity range).

Run directly (not through pytest) so it can drive gz sim as a subprocess
and report a clear pass/fail per case:
    python3 tests/test_continuum_arm_cosserat.py
"""
import math
import os
import subprocess
import sys
import time

from gz.transport13 import Node
from gz.msgs10.vector3d_pb2 import Vector3d
from gz.msgs10.pose_pb2 import Pose

GZ_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORLD = os.path.join(GZ_DIR, "worlds", "continuum_arm_cosserat_test.sdf")

TENDON_TOPIC = "/continuum_arm_cosserat/tendon_cmd"
TIP_POSE_TOPIC = "/continuum_arm_cosserat/tip_pose"
REST_LENGTH = 0.40
PITCH_RADIUS = 0.015
SETTLE_SECONDS = 4.0
POLL_TIMEOUT = 8.0

# l1, l2, l3. +90deg from the textbook 0/120/240deg split: the real
# universal-joint axis convention this arm's plugins share (see
# hillstar_control/continuum_arm_controller.py's estimate_bending_state()
# docstring) bends tendon 0's own commanded shorten toward world +Y, not
# +X -- a joint-rotation-direction quirk verified empirically on
# continuum_arm_pcc, and (after CosseratArmController.cc's OnTendonCmd()
# azimuth fix -- see its own comment) continuum_arm_cosserat now matches
# that same real-world direction for identical tendon commands.
TENDON_ANGLES = [math.pi / 2.0,
                  math.pi / 2.0 + 2.0 * math.pi / 3.0,
                  math.pi / 2.0 + 4.0 * math.pi / 3.0]


def radial_direction(tendon_index):
    theta = TENDON_ANGLES[tendon_index]
    return (math.cos(theta), math.sin(theta))


def read_tip_pose(node, timeout=POLL_TIMEOUT):
    """One-shot subscribe/collect of a single fresh tip_pose message."""
    box = {}

    def cb(msg):
        box["msg"] = msg

    node.subscribe(Pose, TIP_POSE_TOPIC, cb)
    deadline = time.time() + timeout
    while "msg" not in box and time.time() < deadline:
        time.sleep(0.05)
    node.unsubscribe(TIP_POSE_TOPIC)
    if "msg" not in box:
        return None
    p = box["msg"].position
    return (p.x, p.y, p.z)


def publish_tendon_cmd(node, l1, l2, l3):
    pub = node.advertise(TENDON_TOPIC, Vector3d)
    msg = Vector3d()
    msg.x, msg.y, msg.z = l1, l2, l3
    for _ in range(5):
        pub.publish(msg)
        time.sleep(0.05)


def run_direction_case(node, name, tendon_index, shorten_by):
    lengths = [REST_LENGTH, REST_LENGTH, REST_LENGTH]
    lengths[tendon_index] -= shorten_by
    print(f"\n=== case: {name} (l{tendon_index+1} shortened by {shorten_by*1000:.0f}mm) ===")

    publish_tendon_cmd(node, *lengths)
    time.sleep(SETTLE_SECONDS)

    tip = read_tip_pose(node)
    if tip is None:
        print("FAIL: no tip_pose message received")
        return False

    rx, ry = radial_direction(tendon_index)
    lateral_dot = tip[0] * rx + tip[1] * ry
    print(f"measured tip position = ({tip[0]:.5f}, {tip[1]:.5f}, {tip[2]:.5f})  "
          f"radial-projection = {lateral_dot:.5f}")

    # Expect the tip to deflect TOWARD the pulled tendon's own radial
    # direction (lateral_dot > 0) - the shortened side becomes the concave
    # (inside) side of the bend, same as a bimetallic strip curling toward
    # its shorter fiber.
    ok = lateral_dot > 0.001
    print("PASS" if ok else "FAIL: expected deflection toward tendon radial direction")
    return ok


def run_neutral_case(node):
    print("\n=== case: neutral (no tendon deflection) ===")
    publish_tendon_cmd(node, REST_LENGTH, REST_LENGTH, REST_LENGTH)
    time.sleep(SETTLE_SECONDS)
    tip = read_tip_pose(node)
    if tip is None:
        print("FAIL: no tip_pose message received")
        return False
    lateral = math.hypot(tip[0], tip[1])
    print(f"measured tip position = ({tip[0]:.5f}, {tip[1]:.5f}, {tip[2]:.5f})  "
          f"lateral magnitude = {lateral:.5f}")
    ok = lateral < 0.03
    print("PASS" if ok else "FAIL: unexpected lateral deflection with no tendon command")
    return ok


def run_continuity_case(node, tendon_index=0, steps=5, max_shorten=0.03):
    print(f"\n=== case: continuity ramp on tendon {tendon_index+1} "
          f"(0 -> {max_shorten*1000:.0f}mm over {steps} steps) ===")
    magnitudes = []
    rx, ry = radial_direction(tendon_index)

    for step in range(steps + 1):
        shorten = max_shorten * step / steps
        lengths = [REST_LENGTH, REST_LENGTH, REST_LENGTH]
        lengths[tendon_index] -= shorten

        publish_tendon_cmd(node, *lengths)
        time.sleep(SETTLE_SECONDS if step == 0 else 2.0)

        tip = read_tip_pose(node)
        if tip is None:
            print(f"FAIL: no tip_pose message received at step {step}")
            return False

        lateral_dot = tip[0] * rx + tip[1] * ry
        magnitude = lateral_dot  # positive = deflected toward the pulled tendon, as expected
        magnitudes.append(magnitude)
        print(f"  step {step}: shorten={shorten*1000:5.1f}mm  "
              f"tip=({tip[0]:.5f}, {tip[1]:.5f}, {tip[2]:.5f})  "
              f"deflection={magnitude:.5f}m")

    ok = True
    # Monotonic non-decreasing, with a small tolerance for settling noise
    # between adjacent steps (the physical PD-driven joints don't reach the
    # solver's exact target instantaneously).
    for i in range(1, len(magnitudes)):
        if magnitudes[i] < magnitudes[i - 1] - 0.005:
            print(f"FAIL: deflection decreased from step {i-1} to {i} "
                  f"({magnitudes[i-1]:.5f} -> {magnitudes[i]:.5f})")
            ok = False
    if magnitudes[-1] < 0.01:
        print("FAIL: final deflection too small to be a meaningful continuity check")
        ok = False

    print("PASS" if ok else "FAIL")
    return ok


def main():
    env = os.environ.copy()
    env["GZ_SIM_RESOURCE_PATH"] = (
        os.path.join(GZ_DIR, "models") + ":" + os.path.join(GZ_DIR, "worlds")
    )
    env["GZ_SIM_SYSTEM_PLUGIN_PATH"] = ":".join(filter(None, [
        os.path.join(GZ_DIR, "plugins", "cosserat_arm_controller", "build"),
        os.path.join(GZ_DIR, "plugins", "tendon_arm_controller", "build"),
        os.path.join(GZ_DIR, "plugins", "gripper_controller", "build"),
        env.get("GZ_SIM_SYSTEM_PLUGIN_PATH", ""),
    ]))

    print(f"Launching gz sim (headless) with world: {WORLD}")
    server = subprocess.Popen(
        ["gz", "sim", "-s", "-r", WORLD],
        env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )

    try:
        time.sleep(4.0)
        if server.poll() is not None:
            print(server.stdout.read())
            print("FAIL: gz sim server exited early")
            sys.exit(1)

        node = Node()

        results = []
        results.append(("neutral", run_neutral_case(node)))
        results.append(("shorten l1 (bend away from tendon 1)",
                         run_direction_case(node, "shorten l1", 0, 0.03)))
        results.append(("shorten l2 (bend away from tendon 2)",
                         run_direction_case(node, "shorten l2", 1, 0.03)))
        results.append(("shorten l3 (bend away from tendon 3)",
                         run_direction_case(node, "shorten l3", 2, 0.03)))
        results.append(("continuous flexing ramp",
                         run_continuity_case(node)))

        print("\n===== SUMMARY =====")
        all_ok = True
        for name, ok in results:
            print(f"[{'PASS' if ok else 'FAIL'}] {name}")
            all_ok = all_ok and ok

        sys.exit(0 if all_ok else 1)

    finally:
        server.terminate()
        try:
            server.wait(timeout=10)
        except subprocess.TimeoutExpired:
            server.kill()


if __name__ == "__main__":
    main()
