from scipy import linspace
from pylab import plot, show

# Constants
G = 6.674 * 10 ** -11  # gravitational constant
M = 5.974 * 10 ** 24  # mass of earth in kg
m = 7.348 * 10 ** 22  # mass of moon in kg
R = 3.844 * 10 ** 8  # Radius of earth in m
omega = 2.662 * 10 ** -6  # angular velocity of moon
accuracy = 10 ** -8


def f(r):
    return G * M * (R - r) ** 2 - G * m * r ** 2 - omega ** 2 * r ** 3 * (R - r) ** 2


def f_prime(r):
    return -2 * G * M * (R - r) - 2 * G * m * r - 3 * omega ** 2 * r ** 2 * (R - r) ** 2  \
           + 2 * omega ** 2 * r ** 3 * (R - r)


# Plot f to estimate root
# r = linspace(0.75 * R, R, 100)
# f_vals = list(map(f, r))
# plot(r, f_vals, 'o')
# show()

# From the plot, r = 3.2e8 seems like a good starting value

# Solve for the roots of f using Newton's method
def find_root(f, g, start_val, accuracy):
    x = start_val
    delta = f(x) / g(x)
    while abs(delta) > accuracy:
        delta = f(x) / g(x)
        x -= delta

    return x


# Note that f is a rapidly increasing function away from the root
r = find_root(f, f_prime, 3.2 * 10 ** 8, accuracy)
print(r)
print(f(r))
