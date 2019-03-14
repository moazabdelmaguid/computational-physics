from scipy import array, arange, cos, zeros, empty,copy
from pylab import plot, show, xlabel, ylabel
from vpython import sphere, vector, rate

# Constants
m = 1  # mass
k = 6  # spring constant
N = 5  # number of masses
omega = 2
t_0 = 0.0
t_f = 20.0
n = 10000
h = (t_f - t_0) / n

# store the positions and velocities in a vector of length 2N, where the first N elements are the positions,
# and the last N are the corresponding velocities

def f(r, t):
    arr = empty(2 * N, float)
    for i in range(2 * N):
        if i < N:
            arr[i] = r[i + N]
        elif i == N:
            arr[i] =  k / m * (r[1] - r[0]) + 1 / m * cos(omega * t)
        elif i == 2 * N - 1:
            arr[2 * N - 1] = k / m * (r[N - 2] - r[N - 1])
        else:
            arr[i] = k / m * (r[i + 1 - N] - 2 * r[i - N] + r[i - 1 - N])
    return arr


tpoints = arange(t_0, t_f, h)
positions = []
r = zeros(2 * N, float)  # initial displacements and speeds = 0
for t in tpoints:
    positions.append(copy(r[0 : N]))
    k1 = h * f(r, t)
    k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
    k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = h * f(r + k3, t + h)
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

# Plot first mass
# plot(tpoints, array(positions,float)[:, 0])
# show()

# Make animation
spacing = 2
length = N * spacing
masses = empty(N, sphere)
for i in range(N):
    masses[i] = sphere(pos=vector(spacing * i - (length - spacing) / 2, 0, 0), radius = 0.1)

for j in range(len(positions)):
    rate(n//10)
    for i in range(N):
        masses[i].pos = vector(spacing * i - (length - spacing) / 2 + positions[j][i], 0, 0)


