from scipy import array, arange
from pylab import plot, show, xlabel, ylabel

# Constants
sigma = 10
r = 28
b = 8 / 3
t_0 = 0
t_f = 50
x_0 = 0
y_0 = 1
z_0 = 0
N = 100000
h = (t_f - t_0) / N

def f_x(x, y, z):
    return sigma * (y - x)


def f_y(x, y, z):
    return r * x - y - x * z


def f_z(x, y, z):
    return x * y- b * z


def f(r):
    x = r[0]
    y = r[1]
    z = r[2]
    return array([f_x(x, y, z), f_y(x, y, z), f_z(x, y, z)], float)

tpoints = arange(t_0, t_f, h)
xpoints = []
ypoints = []
zpoints = []
R = array([x_0, y_0, z_0], float)
for t in tpoints:
    xpoints.append(R[0])
    ypoints.append(R[1])
    zpoints.append(R[2])
    k1 = h * f(R)
    k2 = h * f(R + 0.5 * h)
    k3 = h * f(R + 0.5 * h)
    k4 = h * f(R + h)
    R += (k1 + 2 * k2 + 2 * k3 + k4)

plot(xpoints, zpoints)
xlabel('t')
ylabel('y(t)')
show()
