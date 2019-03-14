from scipy import array, arange, pi, sin, cos, sqrt
from pylab import plot, show, xlabel, ylabel

# Constants
g = 9.81
# m = 1  # kg
R = 0.08  # m
theta_0 = 30 * pi / 180  # radians
v_0 = 100  # m/s
rho = 1.22  # kg/m^3
C = 0.47  # drag coefficient

t_0 = 0
t_f = 7
N = 10000
h = (t_f - t_0) / N


c = pi * R ** 2 * rho * C / 2
def overall_constant(m):
    return c / m

def f(r, t, m):
    # x = r[0]
    vx = r[1]
    # y = r[2]
    vy = r[3]
    v = sqrt(vx ** 2 + vy ** 2)
    return array([vx, - overall_constant(m) * vx * v,
                  vy, -g - overall_constant(m) * vy * v], float)


tpoints = arange(t_0, t_f, h)
def trajectory(m):
    xpoints = []
    ypoints = []
    r = array([0, v_0 * cos(theta_0), 0, v_0 * sin(theta_0)], float)
    for t in tpoints:
        xpoints.append(r[0])
        ypoints.append(r[2])
        k1 = h * f(r, t, m)
        k2 = h * f(r + 0.5 * k1, t + 0.5 * h, m)
        k3 = h * f(r + 0.5 * k2, t + 0.5 * h, m)
        k4 = h * f(r + k3, t + h, m)
        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return array(xpoints, float), array(ypoints, float)

trajectory1_x, trajectory1_y = trajectory(1)
trajectory2_x, trajectory2_y = trajectory(2)
trajectory3_x, trajectory3_y = trajectory(4)
plot(trajectory1_x, trajectory1_y, 'k')
plot(trajectory2_x, trajectory2_y, 'b')
plot(trajectory3_x, trajectory3_y, 'g')
xlabel('x (m)')
ylabel('y (m)')
show()
