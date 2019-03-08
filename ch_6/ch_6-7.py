from scipy.linalg import solve_banded
from scipy import zeros
from pylab import plot, show

# Constants
N = 10000
V_pos = 5 # volts

# create banded matrix
A = zeros([5, N], float)

A[2, 0] = A[2, N-1] = 3
for i in range(2, N):
    A[0, i] = -1

for i in range(1, N):
    A[1, i] = -1

for i in range(1, N-1):
    A[2, i] = 4

for i in range(N-1):
    A[3, i] = -1

for i in range(N-2):
    A[4, i] = -1

# create RHS vector
b = zeros(N)
b[0] = b[1] = V_pos

# solve system of equations
x = solve_banded((2,2), A, b)

# plot
plot(x, 'ko')
show()
