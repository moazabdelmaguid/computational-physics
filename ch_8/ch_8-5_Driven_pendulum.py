from scipy import sin, cos, array, arange, sqrt
from pylab import plot, show, xlabel, ylabel

# Constants
C = 2  # Hz
l = 0.1  # m
theta_0 = 0.0
omega_0 = 0.0
t_0 = 0.0
t_max = 100  # s
g = 9.81
N = 5000  # number of steps
h = (t_max - t_0) / N

def f(r, t, Omega):
    theta = r[0]
    omega = r[1]
    return array([omega, -(g / l) * sin(theta) + C * cos(theta) * sin(Omega * t)], float)


def theta(Omega):
    tpoints = arange(t_0, t_max, h)
    thetapoints = []
    r = array([theta_0, omega_0], float)
    for t in tpoints:
        thetapoints.append(r[0])
        k1 = h * f(r, t, Omega)
        k2 = h * f(r + 0.5 * k1, t + 0.5 * h, Omega)
        k3 = h * f(r + 0.5 * k2, t + 0.5 * h, Omega)
        k4 = h * f(r + k3, t + h, Omega)
        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return thetapoints


tpoints = arange(t_0, t_max, h)
plot(tpoints, theta(5), 'r')
plot(tpoints, theta(sqrt(g/l)), 'k')
xlabel('t (s)')
ylabel('theta (radians)')
show()
