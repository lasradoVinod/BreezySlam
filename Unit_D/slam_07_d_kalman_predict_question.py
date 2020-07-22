# The complete Kalman prediction step (without the correction step).
#
# slam_07_d_kalman_predict_solution
# Claus Brenner, 12.12.2012
from lego_robot import *
from math import sin, cos, pi, atan2
from numpy import *


class ExtendedKalmanFilter:
    def __init__(self, state, covariance,
                 robot_width,
                 control_motion_factor, control_turn_factor):
        # The state. This is the core data of the Kalman filter.
        self.state = state
        self.covariance = covariance

        # Some constants.
        self.robot_width = robot_width
        self.control_motion_factor = control_motion_factor
        self.control_turn_factor = control_turn_factor

    @staticmethod
    def g(state, control, w):
        x, y, theta = state
        l, r = control
        if r != l:
            alpha = (r - l) / w
            rad = l/alpha
            g1 = x + (rad + w/2.)*(sin(theta+alpha) - sin(theta))
            g2 = y + (rad + w/2.)*(-cos(theta+alpha) + cos(theta))
            g3 = (theta + alpha + pi) % (2*pi) - pi
        else:
            g1 = x + l * cos(theta)
            g2 = y + l * sin(theta)
            g3 = theta

        return array([g1, g2, g3])

    @staticmethod
    def dg_dstate(state, control, w):

        theta = state[2]
        l, r = control
        if r != l:
            alpha = (r-l)/w
            R = l / alpha
            c1 = (R + (w/2)) * (cos (theta + alpha) - cos (theta))
            c2 = (R + (w/2))  * (sin (theta + alpha) - sin (theta))
            m = array([[1, 0, c1], [0, 1, c2], [0, 0, 1]])  # Replace this.
        else:
            m = array([[1, 0, -l * sin(theta)], [0, 1, l * cos(theta)], [0, 0, 1]])  # Replace this.
        return m


    @staticmethod
    def dg_dcontrol(state, control, w):

        theta = state[2]
        l, r = tuple(control)
        dg_3 = -1/w
        dgr_3 = 1/w

        if r != l:
            alpha = (r-l)/w
            R = l / alpha
            first_coeff = (w / ((r - l) ** 2))
            second_coeff = (r+l)/(2 *(r-l))
            dg_1 = first_coeff * r * (sin(theta + alpha) - sin (theta))
            dg_1 = dg_1 - (second_coeff) * cos(theta + alpha)

            dg_2 = first_coeff * r * (cos (theta) - cos(theta + alpha))
            dg_2 = dg_2 - second_coeff * sin (theta + alpha)


            dgr_1 = -first_coeff * l * (sin(theta + alpha) - sin (theta))
            dgr_1 = dgr_1 + second_coeff * cos(theta + alpha)

            dgr_2 = -first_coeff * l * (cos (theta) - cos(theta + alpha))
            dgr_2 = dgr_2 + second_coeff * sin(theta + alpha)

        else:
            dg_1 = 0.5 * (cos (theta) + (l/w) * sin(theta))

            dg_2 = 0.5 * (sin (theta) - (l/w) * cos (theta))

            dgr_1 = 0.5 * (cos (theta) - (l/w) * sin(theta))

            dgr_2 = 0.5 * (sin (theta) + (l/w) * cos(theta))

        m = array([[dg_1, dgr_1], [dg_2, dgr_2], [dg_3, dgr_3]])  # Remove this.
            
        return m
            
        return array([[1, 2], [3, 4], [5, 6]])

    @staticmethod
    def get_error_ellipse(covariance):
        """Return the position covariance (which is the upper 2x2 submatrix)
           as a triple: (main_axis_angle, stddev_1, stddev_2), where
           main_axis_angle is the angle (pointing direction) of the main axis,
           along which the standard deviation is stddev_1, and stddev_2 is the
           standard deviation along the other (orthogonal) axis."""
        eigenvals, eigenvects = linalg.eig(covariance[0:2,0:2])
        angle = atan2(eigenvects[1,0], eigenvects[0,0])
        return (angle, sqrt(eigenvals[0]), sqrt(eigenvals[1]))        

    def predict(self, control):
        """The prediction step of the Kalman filter."""
        # covariance' = G * covariance * GT + R
        # where R = V * (covariance in control space) * VT.
        # Covariance in control space depends on move distance.
        left, right = control

        cov_left = (self.control_motion_factor * left) ** 2
        cov_left = cov_left + (self.control_turn_factor * (left - right)) ** 2


        cov_right = (self.control_motion_factor * right) ** 2
        cov_right = cov_right + (self.control_turn_factor * (left - right)) ** 2

        covariance_control = diag([cov_left,cov_right])
        
        G = self.dg_dstate(self.state,control,self.robot_width)
        V = self.dg_dcontrol(self.state,control,self.robot_width)

        self.covariance = dot(dot(G, self.covariance), transpose(G)) + dot(dot(V, covariance_control), transpose(V))

        # --->>> Put your code to compute the new self.covariance here.

        self.state = self.g(self.state, control, self.robot_width)
    

        # --->>> Put your code to compute the new self.state here.


if __name__ == '__main__':
    # Robot constants.
    scanner_displacement = 30.0
    ticks_to_mm = 0.349
    robot_width = 155.0

    # Filter constants.
    control_motion_factor = 0.35  # Error in motor control.
    control_turn_factor = 0.6  # Additional error due to slip when turning.

    # Measured start position.
    initial_state = array([1850.0, 1897.0, 213.0 / 180.0 * pi])
    # Covariance at start position.
    initial_covariance = diag([100.0**2, 100.0**2, (10.0 / 180.0 * pi) ** 2])
    # Setup filter.
    kf = ExtendedKalmanFilter(initial_state, initial_covariance,
                              robot_width,
                              control_motion_factor, control_turn_factor)

    # Read data.
    logfile = LegoLogfile()
    logfile.read("robot4_motors.txt")

    # Loop over all motor tick records and generate filtered positions and
    # covariances.
    # This is the Kalman filter loop, without the correction step.
    states = []
    covariances = []
    for m in logfile.motor_ticks:
        # Prediction.
        control = array(m) * ticks_to_mm
        kf.predict(control)

        # Log state and covariance.
        states.append(kf.state)
        covariances.append(kf.covariance)

    # Write all states, all state covariances, and matched cylinders to file.
    f = open("kalman_prediction.txt", "w")
    for i in xrange(len(states)):
        # Output the center of the scanner, not the center of the robot.
        print >> f, "F %f %f %f" % \
            tuple(states[i] + [scanner_displacement * cos(states[i][2]),
                               scanner_displacement * sin(states[i][2]),
                               0.0])
        # Convert covariance matrix to angle stddev1 stddev2 stddev-heading form
        e = ExtendedKalmanFilter.get_error_ellipse(covariances[i])
        print >> f, "E %f %f %f %f" % (e + (sqrt(covariances[i][2,2]),))

    f.close()
