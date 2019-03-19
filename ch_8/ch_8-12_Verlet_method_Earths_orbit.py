from scipy import array, arange, sqrt
from pylab import plot, show, xlabel, ylabel

# Constants
x_0 = 1.4710 * 10 ** 11
vx_0 = 0
y_0 = 0
vy_0 = 3.0287 * 10 ** 4 * 8760 * 60 * 60  # m/yr
t_0 = 0
t_f = 4  # year
h = 1 / 8760  # 1 hour in years
G = 6.6738 * 10 ** -11 * ( 8760 * 60 * 60) ** 2
M = 1.9891 * 10 ** 30  # mass of sun in kg
m = 5.9722 * 10 ** 24  # mass of earth in kg


def f(r):
    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]
    dist = sqrt(x ** 2 + y ** 2)
    return array([ vx, -G * M * x / dist ** 3, vy, -G * M * y / dist ** 3 ], float)


# Calculate orbit
tpoints = arange(t_0, t_f, h)
xpoints = []
ypoints = []
potential_E = []
kinetic_E = []
r = array([x_0, vx_0, y_0, vy_0], float)
f_mid = 0.5 * h * f(r)
vx_mid = r[1] + f_mid[1]
vy_mid = r[3] + f_mid[3]
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[2])
    potential_E.append(-6.6738 * 10 ** -11 * M * m / sqrt(r[0] ** 2 + r[2] ** 2))
    kinetic_E.append(0.5 * m * (r[1] ** 2 + r[3] ** 2) / (8760 * 60 * 60) ** 2)
    r[0] += h * vx_mid
    r[2] += h * vy_mid
    k = h * f(r)
    r[1] = vx_mid + 0.5 * k[1]
    r[3] = vy_mid + 0.5 * k[3]
    f_mid = 0.5 * h * f(r)
    vx_mid += k[1]
    vy_mid += k[3]


# Plot orbit
# plot(xpoints, ypoints)
# xlabel('x (m)')
# ylabel('y (m)')
# show()


# Plot energies
total_energy = array(kinetic_E, float) + array(potential_E, float)
plot(tpoints, kinetic_E, 'r')
plot(tpoints, potential_E, 'b')
plot(tpoints, total_energy, 'k')
xlabel('t (years)')
ylabel('Energies (J)')
show()


