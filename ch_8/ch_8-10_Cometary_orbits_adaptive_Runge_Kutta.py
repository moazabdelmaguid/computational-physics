from scipy import array, arange, sqrt, power
from pylab import plot, show, xlabel, ylabel

# Constants
m_sun = 1.989 * 10 ** 30  # kg
G = 66374.2  # in m^3/ kg*yr^2
x_0 = 4 * 10 ** 12  # m
y_0 = 0
v_x = 0
v_y = 15768000000  # m/year
t_0 = 0
t_f = 165  # years

N = 200000
h = (t_f - t_0) / N

def f(r, t):
    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]
    dist = sqrt(x ** 2 + y ** 2)
    return array([vx, - G * m_sun * x / dist ** 3, vy, - G * m_sun * y / dist ** 3], float)


# Using fixed step size
tpoints = arange(t_0, t_f, h)
xpoints = []
ypoints = []
r = array([x_0, v_x, y_0, v_y], float)
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[2])
    k1 = h * f(r, t)
    k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
    k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = h * f(r + k3, t + h)
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

# plot(array(xpoints, float) / 1000, array(ypoints, float) / 1000)
# xlabel('x (km)')
# ylabel('y (km)')
# show()

# Using adaptive step size
def time_step(r, t, h):
    def runge_kutta_step(r, t, h):
        '''
        :param r: current positions and velocities
        :param t: current t
        :param h: step size
        :return: a vector of the change in positions and velocities to get to t+h
        '''
        k1 = h * f(r, t)
        k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
        k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
        k4 = h * f(r + k3, t + h)
        return (k1 + 2 * k2 + 2 * k3 + k4) / 6

    # perform 2 RK steps of step size h
    delta_step_1 = runge_kutta_step(r, t, h)
    delta_step_2 = runge_kutta_step(r + delta_step_1, t + h, h)
    delta_r1 = delta_step_1 + delta_step_2

    # perform 1 RK step with step size 2h
    delta_r2 = runge_kutta_step(r, t, 2 * h)

    # Compute error estimate
    delta_x1 = delta_r1[0]
    delta_x2 = delta_r2[0]
    delta_y1 = delta_r1[2]
    delta_y2 = delta_r2[2]
    error = sqrt((delta_x1 - delta_x2) ** 2 + (delta_y1 - delta_y2) ** 2) / 30

    # Calculate rho
    rho = h * delta / error

    # Calculate factor to multiply h by
    factor = power(rho, 1 / 4)

    # Update h accordingly
    # If target accuracy met, move on to next step
    if  rho >= 1:
        # update t
        t = t + 2 * h

        # Prevent h from getting too large
        if factor > 2:
            h *= 2
        else:
            h *= factor

        # Use local extrapolation to better our estimate of the positions
        delta_r1[0] += (delta_x1 - delta_x2) / 15
        delta_r1[2] += (delta_y1 - delta_y2) / 15
        return delta_r1, h, t
    # If target accuracy not met, must redo step with smaller h
    else:
        return time_step(r, t, factor * h)


delta = 1000   # target accuracy per unit interval in m/yr
h = (t_f - t_0) / 150000  # initial step size
tpoints = []
xpoints2 = []
ypoints2 = []
r = array([x_0, v_x, y_0, v_y], float)  # initial conditions
t = t_0
while(t < t_f):
    tpoints.append(t)
    xpoints2.append(r[0])
    ypoints2.append(r[2])
    delta_r, h, t = time_step(r, t, h)
    r += delta_r


plot(array(xpoints, float) / 1000, array(ypoints, float) / 1000, 'b')
# plot(array(xpoints, float)[::20] / 1000, array(ypoints, float)[::20] / 1000, 'go')
# plot(array(xpoints2, float) / 1000, array(ypoints2, float) / 1000, 'k')
plot(array(xpoints2, float)[::20] / 1000, array(ypoints2[::20], float) / 1000, 'ro')
# plot(tpoints, xpoints, 'ko')
# plot(tpoints, xpoints2, 'g')
# plot(tpoints, ypoints, 'ko')
# plot(tpoints, ypoints2, 'g')
xlabel('x (km)')
ylabel('y (km)')
show()

