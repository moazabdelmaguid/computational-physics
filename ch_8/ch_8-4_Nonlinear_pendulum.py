from scipy import array, arange, pi, sin, cos
from pylab import plot, show, xlabel, ylabel
from vpython import cylinder, vector, sphere, rate

# Constants
g = 9.81
l = 0.1
theta_0 = 179 * pi / 180
omega_0 = 0.0
t_0 = 0.0
t_max = 10.0
N = 5000
h = (t_max - t_0) / N

def f(r, t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = - (g / l) * sin(theta)
    return array([ftheta, fomega], float)


# Using fourth-order Runge-Kutta
tpoints = arange(t_0, t_max, h)
thetapoints = []
# omegapoints = []
r = array([theta_0, omega_0], float)
for t in tpoints:
    thetapoints.append(r[0])
    # omegapoints.append(r[1])
    k1 = h * f(r, t)
    k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
    k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = h * f(r + k3, t + h)
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

# Plot theta vs t
# plot(tpoints, (array(thetapoints, float) * 180 / pi))
# xlabel('t (s)')
# ylabel('theta (degrees)')
# show()

# make animation
rod = cylinder(pos=vector(0, 0, 0), axis=vector(l * cos(theta_0 - pi / 2), l * sin(theta_0 - pi / 2), 0), radius=l/40)
bob = sphere(pos=vector(l * cos(theta_0 - pi / 2), l * sin(theta_0 - pi / 2), 0), radius=l/10)
for theta in thetapoints:
    rate(N // 10)
    rod.axis = vector(l * cos(theta - pi / 2), l * sin(theta - pi / 2), 0)
    bob.pos = vector(l * cos(theta - pi / 2), l * sin(theta - pi / 2), 0)
