from scipy import array, arange, pi, sin, cos
from pylab import plot, show, xlabel, ylabel
from vpython import cylinder, vector, sphere, rate

# Constants
g = 9.81  # m/s^2
m = 1  # kg
l = 0.4  # pendulum lengths in m
theta1_0 = pi / 2
theta2_0 = pi / 2
omega1_0 = 0
omega2_0 = 0
t_0 = 0
t_f = 100
N = 75000
h = (t_f - t_0) / N


# Use fourth-order Runge-Kutta method to solve system
def f(r):
    theta1 = r[0]
    omega1 = r[1]
    theta2 = r[2]
    omega2 = r[3]
    f_omega1 = - (omega1 ** 2 * sin(2 * theta1 - 2 * theta2) + 2 * omega2 ** 2 * sin(theta1 - theta2) + \
                  g / l * (sin(theta1 - 2 * theta2) + 3 * sin(theta1))) / (3 - cos(2 * theta1 - 2 * theta2))
    f_omega2 = (4 * omega1 ** 2 * sin(theta1 - theta2) + omega2 ** 2 * sin(2 * theta1 - 2 * theta2) + \
                 2 * g / l * (sin(2 * theta1 - theta2) - sin(theta2))) / (3 - cos(2 * theta1 - 2 * theta2))
    return array([ omega1, f_omega1, omega2, f_omega2], float)


def energy(r):
    theta1 = r[0]
    omega1 = r[1]
    theta2 = r[2]
    omega2 = r[3]
    return - m * g * l * (2 * cos(theta1) + cos(theta2)) + \
           m * l ** 2 * (omega1 ** 2 + 0.5 * omega2 ** 2 + omega1 * omega2 * cos(theta1 - theta2))


r = array([theta1_0, omega1_0, theta2_0, omega2_0], float)
tpoints = arange(t_0, t_f, h)
theta1_points = []
theta2_points = []
energy_points = []
for t in tpoints:
    theta1_points.append(r[0])
    theta2_points.append(r[2])
    energy_points.append(energy(r))
    k1 = h * f(r)
    k2 = h * f(r + 0.5 * k1)
    k3 = h * f(r + 0.5 * k2)
    k4 = h * f(r + k3)
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

# plot(tpoints, theta1_points,'b')
# plot(tpoints, theta2_points, 'g')
# plot(tpoints, energy_points)
# xlabel('t (s)')
# ylabel('energy (J)')
# show()

# Make animation
rod1 = cylinder(pos=vector(0, 0, 0), axis=vector(l * cos(theta1_0 - pi / 2), l * sin(theta1_0 - pi / 2), 0), radius=l/40)
bob1 = sphere(pos=vector(l * cos(theta1_0 - pi / 2), l * sin(theta1_0 - pi / 2), 0), radius=l/10)
rod2 = cylinder(pos=vector(l * cos(theta1_0 - pi / 2), l * sin(theta1_0 - pi / 2), 0), \
                axis=vector(l * cos(theta2_0 - pi / 2), l * sin(theta2_0 - pi / 2), 0), radius=l/40)
bob2 = sphere(pos=vector(l * cos(theta2_0 - pi / 2), l * sin(theta2_0 - pi / 2), 0), radius=l/10)

for i in range(N):
    rate(N // 100)
    vector1 = vector(l * cos(theta1_points[i] - pi / 2), l * sin(theta1_points[i] - pi / 2), 0)
    vector2 = vector(l * cos(theta2_points[i] - pi / 2), l * sin(theta2_points[i] - pi / 2), 0)
    rod1.axis = vector1
    bob1.pos = vector1
    rod2.pos = vector1
    rod2.axis = vector2
    bob2.pos = vector1 + vector2

