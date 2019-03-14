from scipy import array, arange, sqrt
from pylab import plot, show, xlabel, ylabel

# Constants
G = 1
M = 10
L = 2
x_0 = 1
y_0 = 0.0
vx_0 = 0.0
vy_0 = 1.0
t_0 = 0.0
t_f = 10.0
N = 5000
h = (t_f - t_0) / N

def f(r, t):
    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]
    distance = sqrt(x ** 2 + y ** 2)
    return array([vx, -G * M * x / (distance ** 2 * sqrt(distance ** 2 + L ** 2 / 4)),
                  vy, -G * M * y / (distance ** 2 * sqrt(distance ** 2 + L ** 2 / 4))], float)

tpoints = arange(t_0, t_f, h)
xpoints = []
ypoints = []
r = array([x_0, vx_0, y_0, vy_0], float)
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[2])
    k1 = h * f(r, t)
    k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
    k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = h * f(r + k3, t + h)
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

plot(xpoints, ypoints)
xlabel('x')
ylabel('y')
show()
