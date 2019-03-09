from scipy import zeros, empty, dot, sqrt, array, copy, identity
from scipy import transpose

def eigen_decomp(A, error):
    """
    gives the eigenvectors and eigenvalues of a symmetric matrix A with tolerance error using QR decomposition
    :param A: symmetric matrix
    :param error: positive float
    :return: a list with first element a matrix whose columns are the eigenvectors of A, and the second element
            a diagonal matrix with elements corresponding to the eigenvalues
    """
    N = A.shape[0]

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

    def off_diags_small_enough(a):
        def is_good_enough(x):
            if abs(x) < error:
                return True
            else:
                return False

        for i in range(N):
            for j in range(N):
                if i == j:
                    continue
                else:
                    if not is_good_enough(a[i, j]):
                        return False
        return True

    V = identity(N)
    while(not off_diags_small_enough(A)):
        # Create matrix to store eigenvectors
        Q, R  = QR_decomp(A)

        # update A
        A = dot(R, Q)

        # update V
        V = dot(V, Q)
    return [ V, A ]

# test
A = array([ [1, 4, 8, 4],
            [4, 2, 3, 7],
            [8, 3, 6, 9],
            [4, 7, 9, 2] ], float)
[eigenvectors, eigenvalues] = eigen_decomp(A, 10 ** -6)
print(eigenvectors)
print(eigenvalues)
