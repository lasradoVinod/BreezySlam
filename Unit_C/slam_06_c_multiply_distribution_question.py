# Multiply a distribution by another distribution.
# 06_c_multiply_distribution
# Claus Brenner, 26 NOV 2012
from pylab import plot, show
from distribution import *

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


if __name__ == '__main__':
    arena = (0,1000)

    # Here is our assumed position. Plotted in blue.
    position_value = 400
    position_error = 100
    position = Distribution.gaussian(position_value, position_error)
    plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
         color='b', linestyle='steps')

    # Here is our measurement. Plotted in green.
    # That is what we read from the instrument.
    measured_value = 800
    measurement_error = 100
    measurement = Distribution.gaussian(measured_value, measurement_error)
    plot(measurement.plotlists(*arena)[0], measurement.plotlists(*arena)[1],
         color='g', linestyle='steps')

    # Now, we integrate our sensor measurement. Result is plotted in red.
    position_after_measurement = multiply(position, measurement)

    import operator
    print max(enumerate(position_after_measurement.values), key=operator.itemgetter(1))
    
    plot(position_after_measurement.plotlists(*arena)[0],
         position_after_measurement.plotlists(*arena)[1],
         color='r', linestyle='steps')

    show()
