from scipy import array, arange
from pylab import plot, show, xlabel, ylabel

# Constants
omega = 1
t_0 = 0
t_f = 50
x_0 = 1  # amplitude
v_0 = 0  # initial speed
N = 5000  # number of steps
h = (t_f - t_0) / N

def f_harmonic(r, t):
    x = r[0]
    v = r[1]
    return array([v, - omega ** 2 * x], float)

# harmonic oscillators
tpoints = arange(t_0, t_f, h)
def x_harmonic(amplitude):
    xpoints = []
    r = [amplitude, v_0]
    for t in tpoints:
        xpoints.append(r[0])
        k1 = h * f_harmonic(r, t)
        k2 = h * f_harmonic(r + 0.5 * k1, t + 0.5 * h)
        k3 = h * f_harmonic(r + 0.5 * k2, t + 0.5 * h)
        k4 = h * f_harmonic(r + k3, t + h)
        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return array(xpoints, float)


# anharmonic oscillators
def f_anharmonic(r, t):
    x = r[0]
    v = r[1]
    return array([v, - omega ** 2 * x ** 3], float)

def x_anharmonic(amplitude):
    xpoints = []
    vpoints = []
    r = array([amplitude, v_0], float)
    for t in tpoints:
        xpoints.append(r[0])
        vpoints.append(r[1])
        k1 = h * f_anharmonic(r, t)
        k2 = h * f_anharmonic(r + 0.5 * k1, t + 0.5 * h)
        k3 = h * f_anharmonic(r + 0.5 * k2, t + 0.5 * h)
        k4 = h * f_anharmonic(r + k3, t + h)
        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return array(xpoints, float), array(vpoints, float)

# # Plots of x(t)
# plot(tpoints, x_anharmonic(x_0)[0])
# plot(tpoints, x_anharmonic(2 * x_0)[0])
# xlabel('t (s)')
# ylabel('x (m)')
# show()


# # Phase space plot
# x, v = x_anharmonic(x_0)
# plot(x, v)
# xlabel('x')
# ylabel('v')
# show()

# van der Pol oscillator
# Constants
t_f = 20
N = 10000  # number of steps
h = (t_f - t_0) / N
def g(r, t, mu):
    x = r[0]
    v = r[1]
    return array([v, mu * (1 - x ** 2) * v - omega ** 2 * x], float)

tpoints = arange(t_0, t_f, h)

def x_van_der_pol(mu):
    xpoints = []
    vpoints = []
    r = array([x_0, v_0], float)
    for t in tpoints:
        xpoints.append(r[0])
        vpoints.append(r[1])
        k1 = h * g(r, t, mu)
        k2 = h * g(r + 0.5 * k1, t + 0.5 * h, mu)
        k3 = h * g(r + 0.5 * k2, t + 0.5 * h, mu)
        k4 = h * g(r + k3, t + h, mu)
        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return array(xpoints, float), array(vpoints, float)


mu1_x, mu1_v = x_van_der_pol(1)
mu2_x, mu2_v = x_van_der_pol(2)
mu3_x, mu3_v = x_van_der_pol(4)
plot(mu1_x, mu1_v, 'r')
plot(mu2_x, mu2_v, 'b')
plot(mu3_x, mu3_v, 'g')
xlabel('x')
ylabel('v')
show()