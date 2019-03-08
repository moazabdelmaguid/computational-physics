from scipy import zeros, empty, dot, sqrt, array, copy
from scipy import transpose

def QR_decomp(A):
    """
    returns a list containing arrays Q and R, the QR decomposition of a symmetric matrix A
    :param A: a symmetric array
    :return: a list of arrays
    """
    N = A.shape[0]
    Q = empty([ N, N ], float)
    R = zeros([ N, N ], float)

    # create q vectors
    A_columns = []
    for i in range(N):
        A_columns.append(A[:, i])

    def length_vector(v):
        return sqrt(dot(v, v))

    u_columns = [ copy(A_columns[0]) ]
    q_columns = [ copy(A_columns[0]) / length_vector(A_columns[0]) ]
    for i in range(1, N):
        u_columns.append(copy(A_columns[i]))
        for j in range(i):
            u_columns[i] -= dot(q_columns[j], A_columns[i]) * q_columns[j]
        q_columns.append(u_columns[i] / length_vector(u_columns[i]))

    print(q_columns)
    # Create matrix Q
    for i in range(N):
        Q[:, i] = q_columns[i]

    # Create matrix R
    for i in range(N):
        for j in range(i, N):
            if i == j:
                R[i, i] = length_vector(u_columns[i])
            else:
                R[i, j] = dot(q_columns[i], A_columns[j])

    return [ Q, R ]

# test
A = array([ [1, 4, 8, 4],
            [4, 2, 3, 7],
            [8, 3, 6, 9],
            [4, 7, 9, 2] ], float)
print(QR_decomp(A))