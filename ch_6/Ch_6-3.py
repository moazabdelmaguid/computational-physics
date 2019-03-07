from numpy import array, empty, copy
from pylab import plot, show, xlabel, ylabel, linspace

def LU_decomp(A, v):
    '''
    solves the system of equations associated with the matrix A and vector v
    :param A: square matrix
    :param v: vector
    :return: vector containing the solutions to the system of eqns
    '''
    # define list to store swaps from pivoting
    swaps = []
    
    N = len(v)
    # Gaussian Elimination
    for m in range(N):
        # Partial pivoting
        largest = abs(A[m, m])
        largest_row = m
        for i in range(m + 1, N):
            if abs(A[i, m]) > largest:
                largest = A[i, m]
                largest_row = i
        if largest_row != m:
            current = copy(A[m, :]) # need to use copy because A[m, :] is a reference
            A[m, :] = A[largest_row, :]
            A[largest_row, :] = current

        # Divide by the diagonal element
        div = A[m,m]
        A[m, :] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m + 1, N):
            mult = A[i, m]
            A[i, :] -= mult * A[m, :]
            v[i] -= mult * v[m]

    # Backsubstitution
    x = empty(N, float)
    for m in range(N-1, -1, -1):
        x[m] = v[m]
        for i in range(m+1, N):
            x[m] -= A[m, i] * x[i]

    return x


## Exercise 6.1
A = array([[4, -1, -1, -1],
           [1, -2, 0, 1],
           [-1, 0, 3, -1],
           [1, 1, 1, -4]], float)
v = array([5, 0, 5, 0], float)

print(gaussian_Elim(A, v))

B = array([[0, 1, 4, 1],
           [3, 4, -1, -1],
           [1, -4, 1, 5],
           [2, -2, 1, 3]], float)
w = array([ -4, 3, 9, 7 ], float)
print(gaussian_Elim(B, w))