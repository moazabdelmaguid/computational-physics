from numpy.linalg import solve
from numpy import array, empty


## Exercise 6.4
A = array([[4, -1, -1, -1],
           [1, -2, 0, 1],
           [-1, 0, 3, -1],
           [1, 1, 1, -4]], float)
v = array([5, 0, 5, 0], float)
x = solve(A, v)
print(x)