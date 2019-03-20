from scipy import array, sqrt, power
from pylab import plot, show, xlabel, ylabel

# Constants
G = 1
m1 = 150
m2 = 200
m3 = 250
x1_0 = 3
y1_0 = 1
x2_0 = -1
y2_0 = -2
x3_0 = -1
y3_0 = 1
vx1_0 = 0
vy1_0 = 0
vx2_0 = 0
vy2_0 = 0
vx3_0 = 0
vy3_0 = 0
t_0 = 0
t_f = 10  # arbitrary units


def f(r, t):
    x1 = r[0]
    vx1 = r[1]
    y1 = r[2]
    vy1 = r[3]
    x2 = r[4]
    vx2 = r[5]
    y2 = r[6]
    vy2 = r[7]
    x3 = r[8]
    vx3 = r[9]
    y3 = r[10]
    vy3 = r[11]
    r12 = sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    r23 = sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)
    r13 = sqrt((x1 - x3) ** 2 + (y1 - y3) ** 2)
    return array([vx1, G * m2 * (x2 - x1) / r12 ** 3 + G * m3 * (x3 - x1) / r13 ** 3,
                  vy1, G * m2 * (y2 - y1) / r12 ** 3 + G * m3 * (y3 - y1) / r13 ** 3,
                  vx2, G * m1 * (x1 - x2) / r12 ** 3 + G * m3 * (x3 - x2) / r23 ** 3,
                  vy2, G * m1 * (y1 - y2) / r12 ** 3 + G * m3 * (y3 - y2) / r23 ** 3,
                  vx3, G * m1 * (x1 - x3) / r13 ** 3 + G * m2 * (x2 - x3) / r23 ** 3,
                  vy3, G * m1 * (y1 - y3) / r13 ** 3 + G * m2 * (y2 - y3) / r23 ** 3 ], float)


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

    # Compute error estimates
    # error for star 1
    delta_x1_m1 = delta_r1[0]
    delta_y1_m1 = delta_r1[2]
    delta_x2_m1 = delta_r2[0]
    delta_y2_m1 = delta_r2[2]
    error_m1 = sqrt((delta_x1_m1 - delta_x2_m1) ** 2 + (delta_y1_m1 - delta_y2_m1) ** 2) / 30

    # error for star 2
    delta_x1_m2 = delta_r1[4]
    delta_y1_m2 = delta_r1[6]
    delta_x2_m2 = delta_r2[4]
    delta_y2_m2 = delta_r2[6]
    error_m2 = sqrt((delta_x1_m2 - delta_x2_m2) ** 2 + (delta_y1_m2 - delta_y2_m2) ** 2) / 30

    # error for star 3
    delta_x1_m3 = delta_r1[8]
    delta_y1_m3 = delta_r1[10]
    delta_x2_m3 = delta_r2[8]
    delta_y2_m3 = delta_r2[10]
    error_m3 = sqrt((delta_x1_m3 - delta_x2_m3) ** 2 + (delta_y1_m3 - delta_y2_m3) ** 2) / 30

    # Use the largest error
    error = max(error_m1, error_m2, error_m3)

    # Calculate rhos
    delta = 0.001  # error per unit time
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
            h *= 1.5
        else:
            h *= factor

        # Use local extrapolation to better our estimate of the positions
        delta_r1[0] += (delta_x1_m1 - delta_x2_m1) / 15
        delta_r1[2] += (delta_y1_m1 - delta_y2_m1) / 15
        delta_r1[4] += (delta_x1_m2 - delta_x2_m2) / 15
        delta_r1[6] += (delta_y1_m2 - delta_y2_m2) / 15
        delta_r1[8] += (delta_x1_m3 - delta_x2_m3) / 15
        delta_r1[10] += (delta_y1_m3 - delta_y2_m3) / 15
        return delta_r1, h, t
    # If target accuracy not met, must redo step with smaller h
    else:
        return time_step(r, t, factor * h)


h = (t_f - t_0) / 400000  # initial step size
tpoints = []
xpoints1 = []
ypoints1 = []
xpoints2 = []
ypoints2 = []
xpoints3 = []
ypoints3 = []
r = array([x1_0, vx1_0, y1_0, vy1_0, x2_0, vx2_0, y2_0, vy2_0, x3_0, vx3_0, y3_0, vy3_0], float)  # initial conditions
t = t_0
while(t < t_f):
    tpoints.append(t)
    xpoints1.append(r[0])
    ypoints1.append(r[2])
    xpoints2.append(r[4])
    ypoints2.append(r[6])
    xpoints3.append(r[8])
    ypoints3.append(r[10])
    delta_r, h, t = time_step(r, t, h)
    r += delta_r

plot(xpoints1, ypoints1, 'b')
plot(xpoints2, ypoints2, 'g')
plot(xpoints3, ypoints3, 'r')
xlabel('x')
ylabel('y')
show()
