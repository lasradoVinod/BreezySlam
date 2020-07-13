# Implement the first move model for the Lego robot.
# 02_a_filter_motor
# Claus Brenner, 31 OCT 2012
from math import sin, cos, pi
from pylab import *
from lego_robot import *

def ang_mod(val):
    v = math.fmod(val,2*pi)
    if v < 0:
        return v + (2*pi)
    return v
# This function takes the old (x, y, heading) pose and the motor ticks
# (ticks_left, ticks_right) and returns the new (x, y, heading).
def filter_step(old_pose, motor_ticks, ticks_to_mm, robot_width):
    # Find out if there is a turn at all.
    if motor_ticks[0] == motor_ticks[1]:
        # No turn. Just drive straight.
        x = old_pose[0] + motor_ticks[0] * cos(old_pose[2]) * ticks_to_mm
        y = old_pose[1] + motor_ticks[0] * sin(old_pose[2]) * ticks_to_mm
        # --->>> Implement your code to compute x, y, theta here.
        return (x, y, old_pose[2])

    else:
        # Turn. Compute alpha, R, etc
        alpha = ((motor_ticks[1] - motor_ticks[0]) * ticks_to_mm) / robot_width)

        R = ((motor_ticks[0] * ticks_to_mm/ alpha))
        mul_factor = (R + (robot_width / 2))

        Cx = old_pose[0]- mul_factor * sin(old_pose[2])
        Cy = old_pose[1]+ mul_factor * cos(old_pose[2])
        n_px = Cx + mul_factor * sin(old_pose[2] + alpha)
        n_py = Cy - mul_factor * cos(old_pose[2] + alpha)

        # --->>> Implement your code to compute x, y, theta here.
        return n_px, n_py, ang_mod(old_pose[2] + alpha)


if __name__ == '__main__':
    # Empirically derived conversion from ticks to mm.
    ticks_to_mm = 0.349

    # Measured width of the robot (wheel gauge), in mm.
    robot_width = 150.0

    # Read data.
    logfile = LegoLogfile()
    logfile.read("robot4_motors.txt")

    # Start at origin (0,0), looking along x axis (alpha = 0).
    pose = (0.0, 0.0, 0.0)

    # Loop over all motor tick records generate filtered position list.
    filtered = []
    for ticks in logfile.motor_ticks:
        pose = filter_step(pose, ticks  , ticks_to_mm, robot_width)
        filtered.append(pose)

    # Draw result.
    for pose in filtered:
        print pose
        plot([p[0] for p in filtered], [p[1] for p in filtered], 'bo')
    show()
