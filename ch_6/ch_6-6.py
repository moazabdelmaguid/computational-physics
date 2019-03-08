from vpython import sphere, vector, rate
from scipy import sin, cos, pi, arange, array, zeros, empty
from scipy.linalg import solve_banded
from pylab import plot, show

# Constants
N = 26
C = 1.0
m = 1.0
k = 6.0
omega = 2.0
alpha = 2 * k - m * omega * omega

# construct banded matrix associated with system of equations for input into solve_banded
banded_A = zeros([3, N], float)
for i in range(1, N):
    banded_A[0, i] = -k

for i in range(1, N-1):
    banded_A[1, i] = alpha

for i in range( N-1):
    banded_A[2, i] = -k

banded_A[1,0] = banded_A[1, N-1] = alpha - k

# RHS of equations
b = zeros(N ,float)
b[0] = C

# solve system
x = solve_banded((1,1), banded_A, b)

# make animation

# define position at time t for i-th sphere
def r(A, omega, i, t):
    return 2*i - N + A * cos(omega * t)

s = empty(N, sphere)
for i in range(N):
    s[i] = sphere(pos=vector(r(x[i], omega, i, 0), 0, 0), radius=0.3)

for t in arange(0, 100, 0.1):
    rate(25)
    for i in range(N):
        s[i].pos = vector(r(x[i], omega, i, t), 0, 0)

# s = sphere(pos = vector( 1, 0, 0 ), radius = 0.1)
# for theta in arange(0, 10 * pi, 0.1):
#     rate(30)
#     x = cos(theta)
#     y = sin(theta)
#     s.pos = vector(x, y, 0)
