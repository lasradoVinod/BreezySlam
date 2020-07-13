# Implement the second move model for the Lego robot.
# The difference to the first implementation is:
# - added a scanner displacement
# - added a different start pose (measured in the real world)
# - result is now output to a file, as "F" ("filtered") records.
#
# 02_b_filter_motor_file
# Claus Brenner, 09 NOV 2012
from math import sin, cos, pi
from lego_robot import *
import math
def ang_mod(val):
    v = math.fmod(val,2*pi)
    if v < 0:
        return v + (2*pi)
    return v
# This function takes the old (x, y, heading) pose and the motor ticks
# (ticks_left, ticks_right) and returns the new (x, y, heading).
def filter_step(old_pose, motor_ticks, ticks_to_mm, robot_width,
                scanner_displacement):

    # Find out if there is a turn at all.
    if motor_ticks[0] == motor_ticks[1]:
        # No turn. Just drive straight.
        x = old_pose[0] + motor_ticks[0] * cos(old_pose[2]) * ticks_to_mm
        y = old_pose[1] + motor_ticks[0] * sin(old_pose[2]) * ticks_to_mm
        # --->>> Use your previous implementation.
        # Think about if you need to modify your old code due to the
        # scanner displacement?
        
        return (x, y, old_pose[2])

    else:
        # Turn. Compute alpha, R, etc.
        # Turn. Compute alpha, R, etc
        alpha = ((motor_ticks[1] - motor_ticks[0]) * ticks_to_mm) / robot_width

        R = ((motor_ticks[0] * ticks_to_mm/ alpha))
        mul_factor = (R + (robot_width / 2))

        Cx = old_pose[0] - mul_factor * sin(old_pose[2]) - scanner_displacement * cos(old_pose[2])
        Cy = old_pose[1] + mul_factor * cos(old_pose[2]) - scanner_displacement * sin(old_pose[2])
        n_px = Cx + mul_factor * sin(old_pose[2] + alpha) + scanner_displacement * cos(old_pose[2] + alpha)
        n_py = Cy - mul_factor * cos(old_pose[2] + alpha) + scanner_displacement * sin(old_pose[2] + alpha)


        return n_px, n_py, ang_mod(old_pose[2] + alpha)
        # --->>> Modify your previous implementation.
        # First modify the the old pose to get the center (because the
        #   old pose is the LiDAR's pose, not the robot's center pose).
        # Second, execute your old code, which implements the motion model
        #   for the center of the robot.
        # Third, modify the result to get back the LiDAR pose from
        #   your computed center. This is the value you have to return.


if __name__ == '__main__':
    # Empirically derived distance between scanner and assumed
    # center of robot.
    scanner_displacement = 30.0

    # Empirically derived conversion from ticks to mm.
    ticks_to_mm = 0.349

    # Measured width of the robot (wheel gauge), in mm.
    robot_width = 150.0

    # Measured start position.
    pose = (1850.0, 1897.0, 213.0 / 180.0 * pi)

    # Read data.
    logfile = LegoLogfile()
    logfile.read("robot4_motors.txt")

    # Loop over all motor tick records generate filtered position list.
    filtered = []
    for ticks in logfile.motor_ticks:
        pose = filter_step(pose, ticks, ticks_to_mm, robot_width,
                           scanner_displacement)
        filtered.append(pose)

    # Write all filtered positions to file.
    f = open("poses_from_ticks.txt", "w")
    for pose in filtered:
        print >> f, "F %f %f %f" % pose
    f.close()
