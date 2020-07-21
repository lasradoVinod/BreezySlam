# Instead of moving a distribution, move (and modify) it using a convolution.
# 06_b_convolve_distribution
# Claus Brenner, 26 NOV 2012
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
    return move(Distribution(b.offset,values), a.offset)  # Replace this by your own result.    

if __name__ == '__main__':
    arena = (0,1000)

    # Move 3 times by 20.
    moves = [20] * 50

    # Start with a known position: probability 1.0 at position 10.
    position = Distribution.unit_pulse(0)
    plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
         linestyle='steps')

    # Now move and plot.
    for m in moves:
        move_distribution = Distribution.gaussian(m, 2)
        position = convolve(position, move_distribution)
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             linestyle='steps')

    ylim(0.0, 1.1)
    show()
