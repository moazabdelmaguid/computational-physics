from numpy import zeros
from pylab import imshow, gray, show

# Constants
V = 1.0  # volts
a = 0.1  # grid spacing in centimeters
N = 100  # number of grid points
delta = 10 ** -6  # volts
omega = 0.9

def larger(a, b):
    if a >= b:
        return a
    else:
        return b


# Solve poisson's equation using Gauss-Seidel relaxation method
max_diff = 2 * delta
# initialize array
phi = zeros([N + 1, N + 1], float)
phi[20:81, 20] = V
phi[20:81, 80] = -V
while max_diff > delta:
    #  reset max difference after each complete update of grid
    max_diff = 0.0
    for i in range(N + 1):
        for j in range(N + 1):
            if not i == 0 and not j == 0 and not i == N and not j == N and not (20 <= i <= 80 and j == 20 or j == 80):
                old_phi = phi[i,j]
                new_phi = (1 + omega) * (phi[i + 1, j] + phi[i - 1, j] + phi[i, j + 1] + phi[i, j - 1]) / 4 \
                          - omega * old_phi
                phi[i, j] = new_phi

                # note the largest change in phi in this update of the grid
                max_diff = larger(max_diff, abs(new_phi - old_phi))
    # print("max_diff = ", max_diff)


imshow(phi)
gray()
show()
