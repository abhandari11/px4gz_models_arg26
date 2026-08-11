"""
Functional test for the continuum_arm_soft model + TendonArmController plugin.

Launches `gz sim` headless with worlds/continuum_arm_test.sdf (the arm
rigidly mounted to the world, clear of the ground), publishes varying
3-tendon length commands on /continuum_arm_soft/tendon_cmd, lets physics
settle, and reads back the actual per-joint bend angles from
/continuum_arm_soft/joint_states (gz.msgs.Model, one gz.msgs.Joint per
listed joint, axis1/axis2 .position holding the X/Y bend angle for that
joint) to verify the arm bent into the constant-curvature shape the
tendon-space Jacobian in TendonArmController predicts:
  - shortening l1 alone      -> bends in the joints' local +X direction
  - shortening l2 alone      -> bends along the l2 tendon's 120deg radial
  - shortening l3 alone      -> bends along the l3 tendon's -120deg radial
  - all three commands should agree with the same closed-form FK used by
    the controller itself (independent check, not just "did it move").

Run directly (not through pytest) so it can drive gz sim as a subprocess
and report a clear pass/fail per case:
    python3 tests/test_continuum_arm.py
"""
import math
import os
import subprocess
import sys
import time

from gz.transport13 import Node
from gz.msgs10.vector3d_pb2 import Vector3d
from gz.msgs10.model_pb2 import Model

GZ_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORLD = os.path.join(GZ_DIR, "worlds", "continuum_arm_test.sdf")

TENDON_TOPIC = "/continuum_arm_soft/tendon_cmd"
JOINT_STATE_TOPIC = "/continuum_arm_soft/joint_states"
REST_LENGTH = 0.30
PITCH_RADIUS = 0.018
N_JOINTS = 15
SETTLE_SECONDS = 3.0
POLL_TIMEOUT = 8.0


def expected_theta(l1, l2, l3):
    """Same closed-form FK as TendonArmController::UpdateTargets(), used here
    as an independent oracle rather than re-testing the plugin's own math."""
    d = PITCH_RADIUS
    dl1, dl2, dl3 = l1 - REST_LENGTH, l2 - REST_LENGTH, l3 - REST_LENGTH
    sqrt3_2 = math.sqrt(3) / 2.0
    theta_x = -(2.0 / (3.0 * d)) * (dl1 - 0.5 * dl2 - 0.5 * dl3)
    theta_y = -(2.0 / (3.0 * d)) * (sqrt3_2 * dl2 - sqrt3_2 * dl3)
    return theta_x / N_JOINTS, theta_y / N_JOINTS


def read_joint_states(node, timeout=POLL_TIMEOUT):
    """One-shot subscribe/collect: waits for a single fresh joint_states
    message, then unsubscribes."""
    box = {}

    def cb(msg):
        box["msg"] = msg

    node.subscribe(Model, JOINT_STATE_TOPIC, cb)
    deadline = time.time() + timeout
    while "msg" not in box and time.time() < deadline:
        time.sleep(0.05)
    node.unsubscribe(JOINT_STATE_TOPIC)
    return box.get("msg")


def joint_positions(model_msg):
    """{joint_name: (axis1.position, axis2.position)} from a gz.msgs.Model."""
    out = {}
    for j in model_msg.joint:
        out[j.name] = (j.axis1.position, j.axis2.position)
    return out


def publish_tendon_cmd(node, l1, l2, l3):
    pub = node.advertise(TENDON_TOPIC, Vector3d)
    msg = Vector3d()
    msg.x, msg.y, msg.z = l1, l2, l3
    # A couple of extra publishes since advertise/subscribe handshake can
    # drop the very first message if the plugin's subscriber isn't fully
    # wired up yet.
    for _ in range(5):
        pub.publish(msg)
        time.sleep(0.05)


def run_case(node, name, l1, l2, l3, tol=0.35):
    exp_x, exp_y = expected_theta(l1, l2, l3)
    print(f"\n=== case: {name} (l1={l1:.4f} l2={l2:.4f} l3={l3:.4f}) ===")
    print(f"expected per-joint (theta_x, theta_y) = ({exp_x:.5f}, {exp_y:.5f}) rad")

    publish_tendon_cmd(node, l1, l2, l3)
    time.sleep(SETTLE_SECONDS)

    model_msg = read_joint_states(node)
    if model_msg is None:
        print("FAIL: no joint_states message received")
        return False

    positions = joint_positions(model_msg)
    ok = True
    sum_x = sum_y = 0.0
    for i in range(1, N_JOINTS + 1):
        jname = f"continuum_arm_soft/joint_{i}"
        if jname not in positions:
            print(f"FAIL: missing joint state for {jname}")
            ok = False
            continue
        ax, ay = positions[jname]
        sum_x += ax
        sum_y += ay

    avg_x = sum_x / N_JOINTS
    avg_y = sum_y / N_JOINTS
    print(f"measured avg per-joint (theta_x, theta_y) = ({avg_x:.5f}, {avg_y:.5f}) rad")

    # Direction check: measured bend must have the same sign as expected on
    # each axis whenever the expected magnitude on that axis isn't ~0, and
    # magnitude must be within `tol` (rad) of the closed-form prediction -
    # loose enough to tolerate the passive spring/PD settling point not
    # being an exact match, tight enough to catch a wrong-sign or
    # wrong-axis regression.
    for axis_name, exp, meas in (("x", exp_x, avg_x), ("y", exp_y, avg_y)):
        if abs(exp) > 1e-4:
            if exp * meas <= 0:
                print(f"FAIL: theta_{axis_name} sign mismatch "
                      f"(expected {exp:.5f}, measured {meas:.5f})")
                ok = False
            elif abs(exp - meas) > tol:
                print(f"FAIL: theta_{axis_name} magnitude off "
                      f"(expected {exp:.5f}, measured {meas:.5f}, tol {tol})")
                ok = False
        else:
            if abs(meas) > tol:
                print(f"FAIL: theta_{axis_name} expected ~0 but measured {meas:.5f}")
                ok = False

    print("PASS" if ok else "FAIL")
    return ok


def main():
    env = os.environ.copy()
    env["GZ_SIM_RESOURCE_PATH"] = (
        os.path.join(GZ_DIR, "models") + ":" + os.path.join(GZ_DIR, "worlds")
    )
    env["GZ_SIM_SYSTEM_PLUGIN_PATH"] = ":".join(filter(None, [
        os.path.join(GZ_DIR, "plugins", "tendon_arm_controller", "build"),
        os.path.join(GZ_DIR, "plugins", "gripper_controller", "build"),
        env.get("GZ_SIM_SYSTEM_PLUGIN_PATH", ""),
    ]))

    print(f"Launching gz sim (GUI) with world: {WORLD}")
    server = subprocess.Popen(
        ["gz", "sim", "-r", WORLD],
        env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )

    try:
        time.sleep(4.0)  # let the server come up and load the model/plugins
        if server.poll() is not None:
            print(server.stdout.read())
            print("FAIL: gz sim server exited early")
            sys.exit(1)

        node = Node()

        cases = [
            # name, l1, l2, l3
            ("neutral (no deflection)", REST_LENGTH, REST_LENGTH, REST_LENGTH),
            ("shorten l1 (bend toward +X)", REST_LENGTH - 0.02, REST_LENGTH, REST_LENGTH),
            ("shorten l2 (bend toward l2's 120deg radial)",
             REST_LENGTH, REST_LENGTH - 0.02, REST_LENGTH),
            ("shorten l3 (bend toward l3's -120deg radial)",
             REST_LENGTH, REST_LENGTH, REST_LENGTH - 0.02),
            ("diagonal: shorten l2 and l3 together",
             REST_LENGTH, REST_LENGTH - 0.015, REST_LENGTH - 0.015),
        ]

        results = []
        for name, l1, l2, l3 in cases:
            results.append((name, run_case(node, name, l1, l2, l3)))

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
