from scipy import array, arange
from pylab import plot, show, xlabel, ylabel

# Constants
h = 0.001
x_0 = 1
x_prime_0 = 0
t_0 = 0
t_f = 50

def f(r):
    x = r[0]
    v = r[1]
    return array([ v, v ** 2 - x -5 ], float)

r = array([x_0, x_prime_0] , float)
tpoints = arange(t_0, t_f, h)
xpoints = []
for t in tpoints:
    xpoints.append(r[0])
    r_mid = r + 0.5 * h * f(r)
    r += h * f(r_mid)


plot(tpoints, xpoints)
xlabel('t')
ylabel('x(t)')
show()

