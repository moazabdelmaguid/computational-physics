from scipy import array, sqrt, copy
from pylab import plot, show, xlabel, ylabel

# Constants
G = 6.6738e-11 * ( 8760 * 60 * 60) ** 2
M = 1.9891 * 10 ** 30  # mass of sun in kg
x_0 = 4 * 10 ** 12  # m
vx_0 = 0
y_0 = 0
vy_0 = 500 * 8760 * 60 * 60  # m/yr
delta = 1  # target accuracy of in m/yr
t_0 = 0
t_f = 63.41958396753  # years

def solution(t_0, t_f):
    r = array([x_0, vx_0, y_0, vy_0], float)
    tpoints = [t_0]
    xpoints = [r[0]]
    ypoints = [r[2]]

    def Bulirsch_Stoer_step(r, t, H):
        """
        Computes the nth Richardson extrapolation until the target accuracy is reached
        :param r: vector containing the initial conditions, of the form [x, dx/dt, y, dy/dt]
        :param H: the interval size
        :return: the vector r at t + H
        """

        def modified_midpoint_step(r, n):
            """
            Computes the value of r(t+H) given the initial value of r(t) using the modified midpoint method with n steps
            :param r: array of the form [x, dx/dt, y, dy/dt]
            :param H: interval size
            :param n: number of steps
            :return: array with same form as r
            """

            def f(r):
                """
                Gives the RHS of the system dr/dt = f(r,t)
                :param r: vector of the form (x, dx/dt, y, dy/dt)
                :return: vector of the same form as r
                """

                def fx(x, y):
                    return -G * M * x / sqrt(x ** 2 + y ** 2) ** 3

                def fy(x, y):
                    return -G * M * y / sqrt(x ** 2 + y ** 2) ** 3

                x = r[0]
                vx = r[1]
                y = r[2]
                vy = r[3]
                return array([vx, fx(x, y), vy, fy(x, y)], float)

            r = copy(r)
            h = H / n
            k = r + 0.5 * h * f(r)
            r += h * f(k)
            for i in range(n - 1):
                k += h * f(r)
                r += h * f(k)

            return 0.5 * (r + k + 0.5 * h * f(r))


        def compute_row_n(R1, n):
            """
            Calculates the n-th row (n >= 2) of the Richardson extrapolation table given the (n-1)th row,
             calculates the error on the R_(n, n-1) estimate and compares it to the target accuracy
             and continues to compute the next row until the target accuracy is satisfied
            :param R1: the first richardson estimate R_11, i.e. from the modified midpt method with 1 step
            :param n: the row to compute
            :return: the Richardson estimate of r such that the estimated error is less than the target accuracy
            """

            def R_n_m(m):
                """
                Computes R_n,m
                :param m: integer <= n
                :return: the vector R_n,m
                """
                return R2[m - 2] + (R2[m - 2] - R1[m - 2]) / ((n / (n - 1)) ** (2 * (m - 1)) - 1)

            if n > 8:
                r1 = Bulirsch_Stoer_step(r, t, H / 2)
                return Bulirsch_Stoer_step(r1, t + H / 2, H / 2)
            else:
                # Compute R_n,1
                R2 = [modified_midpoint_step(r, n)]
                # Compute the rest of the row
                for m in range(2, n + 1):
                    R2.append(R_n_m(m))

                # Convert to array to compute error
                R2 = array(R2, float)
                error_vector = (R2[n - 2] - R1[n - 2]) / ((n / (n - 1)) ** (2 * (n - 1)) - 1)
                error = sqrt(error_vector[0] ** 2 + error_vector[2] ** 2)

                # If error is smaller than accuracy, calculation terminates, else repeat with 1 more step
                target_accuracy = H * delta
                if error < target_accuracy:
                    tpoints.append(t + H)
                    xpoints.append(R2[n - 1][0])
                    ypoints.append(R2[n - 1][2])
                    return R2[n - 1]
                else:
                    return compute_row_n(R2, n + 1)


        return compute_row_n(array([modified_midpoint_step(r, 1)], float), 2)

    Bulirsch_Stoer_step(r, t_0, t_f - t_0)
    return tpoints, xpoints, ypoints


t, x, y = solution(t_0, t_f)
plot(x, y, 'o')
xlabel('x (m)')
ylabel('y (m)')
show()

