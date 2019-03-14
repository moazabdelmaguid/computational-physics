from scipy import array, arange
from pylab import plot, show, xlabel, ylabel

# Constants
alpha = 1
beta = 0.5
gamma = 0.5
delta = 2
x_0 = 2
y_0 = 2
t_initial = 0
t_max = 30
N = 20000
h = (t_max - t_initial) / N

def f_x(x, y):
    return alpha * x - beta * x * y


def f_y(x, y):
    return gamma * x * y - delta * y


def f(r):
    x = r[0]
    y = r[1]
    return array([ f_x(x, y), f_y(x, y) ] , float)


tpoints = arange(t_initial, t_max, h)
xpoints = []
ypoints = []

r = array([ x_0, y_0 ], float)
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    k1 = h * f(r)
    k2 = h * f(r + 0.5 * k1)
    k3 = h * f(r + 0.5 * k2)
    k4 = h * f(r + k3)
    r += (k1 + 2 * k2 + 2 * k3 + k4)

plot(tpoints, xpoints, 'b')
plot(tpoints, ypoints, 'r')
xlabel('t')
ylabel('x(t), y(t)')
show()

