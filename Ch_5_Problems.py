from numpy import loadtxt, sum, array, linspace, exp, arange, pi, cos, sin, sqrt, empty, log
from math import factorial, tanh, cosh
from gaussxw import gaussxwab
from pylab import plot, show, xlabel, ylabel, imshow, hot, xlim, ylim, gray

## Exercise 5.1
# data = loadtxt("../cpresources/velocities.txt", float)
# times = data[:, 0]
# velocities = data[:, 1]
#
# def x(t):
#     if t == 0:
#         return 0
#     a = times[0]
#     b = times[int(t)]
#     h = (b - a) / t
#
#     return h*(0.5*velocities[0] + 0.5*velocities[int(t)] + sum(velocities[1:int(t)]))
#
# positions = array(list(map(x, times)),float)
#
# plot(times, velocities)
# plot(times, positions)
# show()

## Exercise 5.2
# a = 0
# b = 2
# N = 1000 # number of slices
# numPoints = N + 1
# # points = linspace(a, b, numPoints)
#
# h = (b-a)/N
#
# def f(x):
#     return x**4 - 2*x + 1
#
# oddSum = 0
# for k in range(1, N, 2):
#     oddSum += f(a + k*h)
#
# evenSum = 0
# for k in range(2, N, 2):
#     evenSum += f(a + k*h)
#
# I = 1/3*h*(f(a) + f(b) + 4*oddSum + 2*evenSum)
# print(I)

# Exercise 5.3
# def f(t):
#     return exp(-t**2)

# integrate using simpson's method
# def E(x):
#     a = 0
#     b = x
#     N = 1000
#     h = (b - a) / N
#     oddSum = 0
#     for k in range(1, N, 2):
#         oddSum += f(a + k*h)
#
#     evenSum = 0
#     for k in range(2, N, 2):
#         evenSum += f(a + k*h)
#
#     return 1 / 3 * h * (f(a) + f(b) + 4 * oddSum + 2 * evenSum)
#
# # make a plot of E(x)
# points = linspace(0, 3, 100)
# values = array(list(map(E, points)))
#
# plot(points, values)
# xlabel("x")
# ylabel("E(x)")
# show()

# # Exercise 5.4
# def J(m, x):
#     def f(m, x, theta):
#         return cos(m*theta - x* sin(theta))
#
#     N = 1000
#     a = 0
#     b = pi
#     h = (b - a) / N
#
#     oddSum = 0
#     for k in range(1, N, 2):
#         oddSum += f(m, x, a + k*h)
#
#     evenSum = 0
#     for k in range(1, N, 2):
#         evenSum += f(m, x, a + k*h)
#
#     return 1 / pi * 1 / 3 * h * (f(m, x, a) + f(m, x, b) + 4 * oddSum + 2 * evenSum)
#
# # Plot J0, J1, J2
# # xpoints = linspace(0, 20, 100)
# # J0 = []
# # J1 = []
# # J2 = []
# # for x in xpoints:
# #     J0.append(J(0, x))
# #     J1.append(J(1, x))
# #     J2.append(J(2, x))
#
# # plot(xpoints, J0, "g")
# # plot(xpoints, J1, "b")
# # plot(xpoints, J2, "r")
# # show()
#
# # Part b
# def r(x, y):
#     return sqrt(x**2 + y**2)
#
# def I(r):
#     if (r == 0):
#         return 1/4
#
#     Lambda = 0.5 # in micrometers
#     kr = 2 * pi / Lambda * r
#     return (J(1,kr)/ kr)**2
#
# side = 2 # length in micrometers
# points = 200 # number of grid points in each direction
# spacing = side/points
#
# # Calculate the position of the center
# xCenter = side/2
# yCenter = side/2
#
# # Make an empty array to store values
# intensities = empty([points, points], float)
#
# # Calculate the values in the array
# for i in range(points):
#     y = spacing * i
#     for j in range(points):
#         x = spacing * j
#         dist = r(x - xCenter, y - yCenter)
#         intensities[i, j] = I(dist)
#
# imshow(intensities, origin="lower", extent=[0,side,0,side], vmax=0.01)
# hot()
# show()

# Exercise 5.6
# def f(x):
#     return x**4 - 2*x + 1
#
# N1 = 10
# N2 = 2*N1
# numPoints1 = N1 + 1
# numPoints2 = N2 + 1
# a = 0
# b = 2
# h1 = (b - a) / N1
# h2 = (b - a) / N2
#
# fValues1 = sum(list(map(f, arange(a + h1, b, h1))))
# fValues2 = sum(list(map(f, arange(a + h2, b, h1)))) # we can nest the points
#
# I1 = h1 * (0.5*f(a) + 0.5*f(b) + fValues1)
# I2 = h2 * (0.5*f(a) + 0.5*f(b) + fValues1 + fValues2)
# epsilon = 1/3 * (I2-I1)
# print(I1)
# print(I2)
# print("error is approximately", epsilon)
# print(I2 + epsilon) # error is good to only h^2

# Exercise 5.7 part a
# def adapTrapI(f, a, b, N, Iprev, error):
#     N = 2*N
#     h = (b - a)/N
#     I = 0.5*Iprev + h*sum(list(map(f, arange(a + h, b, 2*h))))
#     print(I)
#     epsilon = 1/3*(I - Iprev)
#     if (abs(epsilon) < error):
#         print('N = ', N)
#         print('error = ', epsilon)
#         return I
#     else:
#         return adapTrapI(f, a, b, N, I, error)
#
# def f(x):
#     return sin(sqrt(100*x))**2
#
# print(adapTrapI(f, 0, 1, 1, 0.5*(f(1)-f(0)), .000001))

# part b
#
# def adapTrapI_step(f, a, b, N, Iprev):
#     N = 2 * N
#     h = (b - a)/N
#     I = 0.5 * Iprev + h * sum(list(map(f, arange(a + h, b, 2 * h))))
#     return I
#
# def romberg(f, a, b, error):
#     #starts with N=2 slices
#     R11 = (b-a)/2 * (0.5 * f(a) + 0.5 * f(b) + f(a + (b-a)/2 )) # initial estimate of I with N=2
#     R_list = [ R11, adapTrapI_step(f, a, b, 2, R11) ] # keep a list of the romberg estimates
#
#     def Rij(i, m):
#         # returns the linear index of Rij corresponding to i, m
#         return R_list[int(1/2 * i * (i - 1) + m - 1)] # given by sum_{j=0}^{i-1} + m - 1
#
#     def romberg_step(i, m):
#         #returns R_i,m such that error estimate is less than given error
#         #print(R_list)
#         if m == 1:
#             # if m = 1, then need to compute the next nested trapezoidal rule estimate
#             print('i = ', i)
#             print('m = ', m)
#             R_list.append(adapTrapI_step(f, a, b, 2**(i-1), Rij(i - 1, 1)))
#             epsilon = 1/3 * (Rij(i, 1) - Rij(i - 1, 1))
#             print(Rij(i, m))
#             print('epsilon = ', epsilon)
#             if abs(epsilon) < error:
#                 return Rij(i, m)
#             else:
#                 return romberg_step(i, m + 1)
#         else:
#             print('i = ', i)
#             print('m = ', m)
#             epsilon = 1/(4**(m-1) - 1) * (Rij(i, m - 1) - Rij(i - 1, m - 1))
#             print('epsilon = ', epsilon)
#             R_list.append(Rij(i, m - 1) + epsilon)
#             print(Rij(i,m))
#             if abs(epsilon) < error:
#                 return Rij(i, m)
#             else:
#                 if i == m:
#                     return romberg_step(i + 1, 1)
#                 else:
#                     return romberg_step(i, m + 1)
#     return romberg_step(2, 2)
#
# def f(x):
#     return sin(sqrt(100*x))**2
#
# def g(x):
#     return x**2
#
# print(romberg(f, 0, 1, 0.000001))


## Exercise 5.9
# def cV(T):
#     c = 7.48279 # = 9V*rho*k_B in SI units
#     thetaD = 428 # in K
#
#     def f(x):
#         return x**4 * exp(x) / (exp(x) - 1)**2
#
#     # perform the integration using Gaussian quadrature
#     N = 50  # 50 sample points
#     x, w = gaussxwab(N, 0, thetaD / T)
#     integral = 0.0
#     for k in range(N):
#         integral += w[k] * f(x[k])
#
#     return c * (T / thetaD) ** 3 * integral
#
# T = linspace(5, 500, 99)
# C = list(map(cV, T))
# plot(T, C, 'o')
# xlim(5, 500)
# xlabel('T (K)')
# ylabel('C_V (J/K)')
# show()


## Exercise 5.10
# m = 1 # mass
# N = 20 # number of points for gaussian quadrature
# def T(a):
#     def f(x):
#         def v(y):
#             return y ** 4
#
#         return 1 / sqrt(v(a) - v(x))
#
#     x, w = gaussxwab(N, 0, a)
#     integral = 0.0
#     for k in range(N):
#         integral += w[k] * f(x[k])
#
#     return sqrt(8)*integral
#
# # make a plot of T for a from 0 to 2
# a = linspace(0, 2, 20)
# periods = list(map(T, a))
# plot(a, periods, 'o')
# xlabel('a (m)')
# ylabel('T (s)')
# show()

# ## Exercise 5.11
# def integrator(f, a, b, N):
#     """
#     integrates f from a to b using Gaussian quadrature with N points, requires gaussxwab
#     :param f: 1d function
#     :param a: lower bound of domain
#     :param b: upper bound of domain
#     :param N: number of points
#     :return:  floating value of integral
#     """
#     x, w = gaussxwab(N, a, b)
#     integral = 0.0
#     for k in range(N):
#         integral += w[k] * f(x[k])
#
#     return integral
#
# def I(x, wavelength, z):
#     '''
#
#     :param x: x in meters
#     :param wavelength: in meters
#     :param z: in meters
#     :return: value of I / I0
#     '''
#     u = x * sqrt(2 / (wavelength * z))
#     N = 50 # number of integration points
#     def c(u):
#         def f(x):
#             return cos(0.5 * pi * x ** 2)
#         return integrator(f, 0, u, N)
#
#     def s(u):
#         def g(x):
#             return sin(0.5 * pi * x ** 2)
#         return integrator(g, 0, u, N)
#
#     return 1 / 8 * ((2 * c(u) + 1) ** 2 + (2 * s(u) + 1) ** 2)
# # plot I/I0 for x between -5 and 5 m, wavelength = 1, z = 3 m
# xvalues = linspace(-5, 5, 100)
# yvalues = []
# wavelength = 1
# z = 3
# for k in range(100):
#     yvalues.append(I(xvalues[k], wavelength, z))
#
# plot(xvalues, yvalues)
# xlabel('x (m)')
# ylabel('I/I0')
# show()


## Exercise 5.12
# def f(z):
#     return (z / (1 - z) ) ** 3 / (exp(z / (1 - z)) - 1) * 1 / (1 - z) ** 2
# # integrand is pretty smooth, so can use gaussian quadrature
# N = 50 # points for Gaussian quadrature
# x, w = gaussxwab(N, 0, 1)
# integral = 0.0
# for k in range(N):
#     integral += w[k] * f(x[k])
# print(integral)

# ## Exercise 5.13
# def H(n, x):
#     '''
#     computes the n-th Hermite polynomial using a linear iterative procedure
#     :param n: integer
#     :param x: value of x
#     :return: value of H_n(x)
#     '''
#     def H_iter(a, b, count):
#         if (count == 0):
#             return b
#         else:
#             return H_iter(2  * x * a - 2 * (count - 1) * b, a, count - 1)
#
#     return H_iter(2 * x, 1, n)
#
# # plot H for n = 0 - 3, x between -4, 4
# # xvalues = linspace(-4, 4, 100)
# # H0values = []
# # H1values = []
# # H2values = []
# # H3values = []
# # for k in range(100):
# #     H0values.append(H(0,xvalues[k]))
# #     H1values.append(H(1, xvalues[k]))
# #     H2values.append(H(2, xvalues[k]))
# #     H3values.append(H(3, xvalues[k]))
# #
# # plot(xvalues, H0values)
# # plot(xvalues, H1values)
# # plot(xvalues, H2values)
# # plot(xvalues, H3values)
# # show()
#
# def psi(n, x):
#     return 1 / sqrt(2 ** n * factorial(n) * sqrt(pi)) * exp(- x ** 2 / 2) * H(n, x)
#
# # xvalues = linspace(-10, 10, 500)
# # H30values = []
# # for k in range(500):
# #     H30values.append(psi(30, xvalues[k]))
# # plot(xvalues, H30values)
# # show()
#
# def rms_integrand(z):
#     def x(z):
#         return z / (1 - z)
#
#     return x(z) ** 2 * abs(psi(5,x(z))) ** 2 * (1/(1 - z) ** 2)
#
# integral = 0.0
# N = 100
# x, w = gaussxwab(N, 0,1)
# for k in range(N):
#     integral += w[k] * rms_integrand(x[k])
# print(sqrt(2*integral))

# ## Exercise 5.14
# G = 6.674 * 10 ** -11 # gravitational constant in SI units
# sigma = 100 # mass density in SI units
# N = 100 # number of points in each direction for double Gaussian quadrature
# r, w = gaussxwab(N, -5,5)
# def Fz(z):
#     def integrand(x, y, z):
#         return 1 / (x**2 + y ** 2 + z ** 2) ** 3/2
#
#     integral = 0
#     for i in range(N):
#         for j in range(N):
#             integral += w[i] * w[j] * z * integrand(r[i], r[j], z)
#
#     return G * sigma * integral
#
# zvalues = linspace(0, 10, 100)
# Fzvalues = list(map(Fz, zvalues))
#
# plot(zvalues, Fzvalues, 'o')
# xlabel('z (m)')
# ylabel('Fz (N)')
# show()
# # The force seems to drop to zero for small values of z because the integrand becomes very small for points
# # far from x,y = 0, and we don't have enough points near (0,0) when calculating the integral

# ## Exercise 5.15
# def f(x):
#     return 1 + 0.5 * tanh(2*x)
#
# # calculate df/dx using central difference method
# def df_dx(x):
#     h = 10 ** -5  # step size
#     return (f(x + 0.5 * h) - f(x - 0.5 * h)) / h
#
# def g(x):
#     # analytic derivative of f(x) above
#     return 1 / cosh(2*x) ** 2
#
# xvals = linspace(-2, 2, 100)
# dfvals = list(map(df_dx, xvals))
# gvals = list(map(g, xvals))
#
# plot(xvals, dfvals, 'o')
# plot(xvals, gvals)
# xlabel('x')
# show()


# ## Exercise 5.17
# # for change of variables z = x / (c + x), x = c gives z = 1/2
# # thus since the max of x^(a-1) e^-x occurs at a-1, choosing c = a-1 puts the peak of the integrand at z = 1/2
# # let's use gaussian quadrature with 100 points
# N = 100
# x, w = gaussxwab(N, 0 ,1)
# def gamma(a):
#     c = a - 1
#     def integrand(z):
#         return c * exp(c * log((c * z) / (1 - z)) - (c * z) / (1 - z)) / (1 - z) ** 2
#     integral = 0
#     for k in range(N):
#         integral += w[k] * integrand(x[k])
#
#     return integral
#
# print(gamma(3/2)) # exact value is sqrt(pi)/2
# print(gamma(10)) # equals 9! = 362880


# ## Exercise 5.19
# focal_len = 1 # in m
# screen_width = 0.1 # in m
# wavelength = 0.5 # in um (micrometers)
# slit_sep = 20 # in um
# num_slits = 10
# grating_width = slit_sep * num_slits
# def q(u):
#     alpha = pi / slit_sep
#     return sin(alpha * u) ** 2
#
# def I(x):
#     # The integrand is highly oscillatory, so let's just use the trapezoidal rule
#     # note we've scaled u in terms of um
#     def integrand(u):
#         return sqrt(q(u)) * exp(2j * pi * x * u / (wavelength * focal_len))
#
#     N = 1000 # num of integration slices
#     h = grating_width / N # step size
#     integral = h * 0.5 *(integrand(- grating_width / 2) + integrand(grating_width / 2))
#     for k in range(1, N, 2):
#         integral += integrand(-grating_width / 2 + k * h)
#
#     return 10 ** -12 * abs(integral) ** 2
#
# xvals = linspace(-.05, .05, 500)
# Ivals = list(map(I, xvals))
#
# # plot(xvals, Ivals, 'o')
# # # show()
#
# Iarray = empty([100,500], float)
# for k in range(100):
#     Iarray[k, :] = Ivals
# imshow(Iarray)
# gray()
# show()


## Exercise 5.20
def f(x):
    if x == 0:
        return 1
    else:
        return sin(x) ** 2 / x ** 2

def integral(f, a, b, error):
    delta = error / (b - a) # target accuracy per unit interval
    def step(x1, x2, f1, f2):
        # calculates estimates of the integral from x1 to x2 with one and two slices, and the est. error
        h = x2 - x1
        midpoint = 0.5 * (x2 + x1)
        f_mid = f(midpoint)
        I1 = h * 0.5 * (f1 + f2)
        I2 = 0.5 * I1 + 0.5 * h * f_mid
        if (abs(1 / 3 * (I2 - I1)) < h * delta):
            return 1 / 6 * h * (f1 + 4 * f_mid + f2)
        else:
            return step(x1, midpoint, f1, f_mid) + step(midpoint, x2, f_mid, f2)

    return step(a, b, f(a), f(b))

print(integral(f, 0, 10, 10 ** -4))




