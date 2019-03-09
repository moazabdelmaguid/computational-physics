from scipy import pi, empty, dot, sqrt, linspace, array, sin
from scipy.linalg import eigh, eigvalsh
from pylab import plot, show, xlabel, ylabel

# Constants
hbar = 6.582119514 * 10 ** -16 # in J * s
L = 5 * 10 ** -10 # in m
M = 9.1094 * 10 ** -31 # in kg
a = 10 # in eV


def h(m, n):
    def is_same_parity(m, n):
        return m % 2 == n % 2

    if m == n:
        return (pi ** 2 * hbar ** 2 * n ** 2) / (2 * M * L ** 2) * (1.6022 * 10 ** -19) + a / 2
    elif not is_same_parity(m, n):
        return - 8 * a / (pi ** 2)  * m * n / (m ** 2 - n ** 2) ** 2
    else:
        return 0


# Create H up to 10 x 10 elements
# H = empty([ 10, 10 ])
# for m in range(10):
#     for n in range(10):
#         H[m, n] = h(m + 1, n + 1)
#
# E, psi = eigh(H)
# print(E)

# Compare with 100 x 100
H = empty([ 100, 100 ], float)
for m in range(100):
    for n in range(100):
        H[m, n] = h(m + 1, n + 1)

E, psi = eigh(H)

# plot three lowest energy states

def V(x):
    return a * x / L

def psi_n(n, x):
    psi0 = 0
    for m in range(100):
        psi0 += sqrt(2 / L) * psi[n][m] * sin(pi * (m+1) * x / L)
    return psi0

def psi_0(x):
    return psi_n(0, x)

def psi_1(x):
    return psi_n(1, x)

def psi_2(x):
    return psi_n(2, x)

def square(x):
    return x ** 2

x = linspace(0, L, 100)
v = list(map(V, x))
psi0 = array(list(map(psi_0, x)))
psi1 = array(list(map(psi_1, x)))
psi2 = array(list(map(psi_2, x)))
psi0_squared = list(map(square, psi0))
psi1_squared = list(map(square, psi1))
psi2_squared = list(map(square, psi2))
plot(x, psi0_squared, 'ko')
plot(x, psi1_squared, 'go')
plot(x, psi2_squared, 'bo')
plot(x, v)
xlabel('x (m)')
show()


