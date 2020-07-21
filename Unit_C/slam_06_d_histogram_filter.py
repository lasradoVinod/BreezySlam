# Histogram implementation of a bayes filter - combines
# convolution and multiplication of distributions, for the
# movement and measurement steps.
# 06_d_histogram_filter
# Claus Brenner, 28 NOV 2012
from pylab import plot, show, ylim
from distribution import *

def move(distribution, delta):
    """Returns a Distribution that has been moved (x-axis) by the amount of
       delta."""
    return Distribution(distribution.offset + delta, distribution.values)

def convolve(a, b):
    length = len(a.values) + len(b.values) - 1
    values = []
    for i in range(0,length):
        value = 0
        for j in range(len(a.values)):
            k = i-len(a.values) + 1 + j
            if k < 0:
                continue
            if k >= len(b.values):
                break      
            value += (a.values[j] * b.values[k])
        values.append(value)
        print a.offset
        print b.offset
    return move(Distribution(b.offset,values), a.offset)  # Replace this by your own result.


def multiply(a, b):
    """Multiply two distributions and return the resulting distribution."""
    res_offset = max(a.offset,b.offset)
    if (res_offset == b.offset):
        start = b
        other = a
    else:
        start = a
        other = b

    values = []
    start_offset = res_offset - other.offset 
    for i in range(len(start.values)):
        if (i + start_offset >= len(other.values)):
            break 
        values.append(start.values[i] * other.values[start_offset + i])

    sum_val = sum(values)
    norm_factor = 1
    if not sum_val == 0:
        norm_factor = 1/sum(values)
    values = [entry * norm_factor for entry in values]
    return Distribution(res_offset,values) # Modify this to return your result.


# --->>> Copy your convolve(a, b) and multiply(a, b) functions here.



if __name__ == '__main__':
    arena = (0,220)

    # Start position. Exactly known - a unit pulse.
    start_position = 10
    position = Distribution.unit_pulse(start_position)
    plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
         linestyle='steps')

    # Movement data.
    controls  =    [ 20 ] * 10

    # Measurement data. Assume (for now) that the measurement data
    # is correct. - This code just builds a cumulative list of the controls,
    # plus the start position.
    p = start_position
    measurements = []
    for c in controls:
        p += c
        measurements.append(p)

    # This is the filter loop.
    for i in xrange(len(controls)):
        # Move, by convolution. Also termed "prediction".
        control = Distribution.triangle(controls[i], 20)
        position = convolve(position, control)
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             color='b', linestyle='steps')

        # Measure, by multiplication. Also termed "correction".
        measurement = Distribution.triangle(measurements[i], 20)
        position = multiply(position, measurement)
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             color='r', linestyle='steps')

    show()
