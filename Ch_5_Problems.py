from numpy import loadtxt, sum, array, linspace, exp, arange, pi, cos, sin, sqrt, empty
from pylab import plot, show, xlabel, ylabel, imshow, hot

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
def f(x):
    return x**4 - 2*x + 1

N1 = 10
N2 = 2*N1
numPoints1 = N1 + 1
numPoints2 = N2 + 1
a = 0
b = 2
h1 = (b - a) / N1
h2 = (b - a) / N2

fValues1 = sum(list(map(f, arange(a + h1, b, h1))))
fValues2 = sum(list(map(f, arange(a + h2, b, h1)))) # we can nest the points

I1 = h1 * (0.5*f(a) + 0.5*f(b) + fValues1)
I2 = h2 * (0.5*f(a) + 0.5*f(b) + fValues1 + fValues2)
epsilon = 1/3 * (I2-I1)
print(I1)
print(I2)
print("error is approximately", epsilon)
print(I2 + epsilon) # error is good to only h^2

