# This adds the derivative of the measurement function h,
# with respect to the state.
#
# slam_07_e_measurement_derivative
# Claus Brenner, 12.12.2012
from lego_robot import *
from math import sqrt, atan2, sin, cos
from numpy import *


class ExtendedKalmanFilter:

    @staticmethod
    def h(state, landmark, scanner_displacement):
        """Takes a (x, y, theta) state and a (x, y) landmark, and returns the
           measurement (range, bearing)."""
        dx = landmark[0] - (state[0] + scanner_displacement * cos(state[2]))
        dy = landmark[1] - (state[1] + scanner_displacement * sin(state[2]))
        r = sqrt(dx * dx + dy * dy)
        alpha = (atan2(dy, dx) - state[2] + pi) % (2*pi) - pi

        return array([r, alpha])

    @staticmethod
    def dh_dstate(state, landmark, scanner_displacement):
        x, y, theta = state
        x_l = x + scanner_displacement * cos(theta)
        y_l = y + scanner_displacement * sin(theta)

        del_x = landmark[0] - x_l
        del_y = landmark[1] - y_l

        q = del_x ** 2 + del_y ** 2
        sqrt_q = sqrt(q)

        drx = -del_x / sqrt_q
        dry = -del_y / sqrt_q
        drt = (scanner_displacement/sqrt_q) * (del_x * sin(theta) - del_y * cos(theta))

        dax = del_y /q
        day = -del_x /q
        dat = -1 + (-scanner_displacement/q) * (del_x * cos(theta) + del_y * sin(theta)) 
     
        return array([[drx, dry, drt], [dax, day, dat]]) # Replace this.


if __name__ == '__main__':
    # If the partial derivative with respect to x, y, theta (the state)
    # are correct, then the numerical derivative and the analytical
    # derivative should be the same.

    # Set some variables. Try other variables as well.
    x = 11.0
    y = 22.0
    theta = 125. / 180. * pi
    state = array([x, y, theta])
    landmark = (203.0, 20.0)
    scanner_displacement = 30.0
    w = 150.0

    # Compute derivative numerically.
    print "Numeric differentiation dx, dy, dtheta:"
    delta = 1e-7
    state_x = array([x + delta, y, theta])
    state_y = array([x, y + delta, theta])
    state_theta = array([x, y, theta + delta])
    dh_dx = (ExtendedKalmanFilter.h(state_x, landmark, scanner_displacement) -\
             ExtendedKalmanFilter.h(state, landmark, scanner_displacement)) / delta
    dh_dy = (ExtendedKalmanFilter.h(state_y, landmark, scanner_displacement) -\
             ExtendedKalmanFilter.h(state, landmark, scanner_displacement)) / delta
    dh_dtheta = (ExtendedKalmanFilter.h(state_theta, landmark, scanner_displacement) -\
                 ExtendedKalmanFilter.h(state, landmark, scanner_displacement)) / delta
    dh_dstate_numeric = column_stack([dh_dx, dh_dy, dh_dtheta])
    print dh_dstate_numeric

    # Use the above code to compute the derivative analytically.
    print "Analytic differentiation dx, dy, dtheta:"
    dh_dstate_analytic = ExtendedKalmanFilter.dh_dstate(
        state, landmark, scanner_displacement)
    print dh_dstate_analytic

    # The difference should be close to zero (depending on the setting of
    # delta, above).
    print "Difference:"
    print dh_dstate_numeric - dh_dstate_analytic
    print "Seems correct:", allclose(dh_dstate_numeric, dh_dstate_analytic)
