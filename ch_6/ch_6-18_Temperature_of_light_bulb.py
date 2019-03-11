from scipy import sqrt, pi, exp, linspace
from gaussxw import gaussxwab
from pylab import plot, show, xlabel, ylabel

# Constants
z = (1 + sqrt(5)) / 2  # golden ratio
hc = 1.23984193 * 10 ** 3  # in eV * nm
lambda_1 = 390  # in nm
lambda_2 = 750  # in nm
kB = 8.6173303 * 10 ** -5  # in eV / K
lower_constant = hc / (lambda_2 * kB)
upper_constant = hc / (lambda_1 * kB)
accuracy = 1  # in K

def eta(T):
    def f(x):
        return x ** 3 / (exp(x) - 1)

    # We use Gaussian quadrature with N = 100 sample points
    N = 100
    x, w = gaussxwab(N, lower_constant / T, upper_constant / T)
    integral = 0.0
    for k in range(N):
        integral += w[k] * f(x[k])

    return 15 / pi ** 4 * integral

# Make plot of eta from T = 300 K to 10000 K
# T = linspace(300, 10000, 100)
# etas = list(map(eta, T))
# plot(T, etas)
# xlabel('T (K)')
# ylabel(('eta'))
# show()

# From plot, max efficiency near T = 7000 K
# Use golden ratio search to calculate temperature of maximum efficiency
# Initial points
T1 = 6000
T4 = 8000
T2 = T4 - (T4 - T1) / z
T3 = T4 + (T4 - T1) / z

# Initial values of eta
eta_1 = eta(T1)
eta_2 = eta(T2)
eta_3 = eta(T3)
eta_4 = eta(T4)

# golden ratio search loop
while T4 - T1 > accuracy:
    if eta_2 < eta_3 :
        T4, eta_4 = T3, eta_3
        T3, eta_3 = T2, eta_2
        T2 = T4 - (T4 - T1) / z
        eta_2 = eta(T2)
    else:
        T1, eta_1 = T2, eta_2
        T2, eta_2 = T3, eta_3
        T3 = T1 + (T4 - T1) / z
        eta_3 = eta(T3)

print('The temperature of max efficiency is', 0.5 * (T1 + T4))

# Considering the melting point of tungsten is 3695 K, you can't run a tungsten-filament light bulb at this temp
