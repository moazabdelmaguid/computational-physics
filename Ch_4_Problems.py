from pylab import plot, show, ylim, xlabel, ylabel, imshow, gray, colorbar, scatter, hot
from numpy import linspace, sin, cos, loadtxt, pi, exp, sqrt, empty, arange, zeros, array, meshgrid, full, dot
import math


# Exercise 4.1
# def factorialInt(n):
#     prod = 1
#     for k in range(n, 0, -1):
#         prod *= k
#     return prod
#
#
# def factorialFloat(n):
#     prod = 1.0
#     for k in range(n, 0, -1):
#         prod *= k
#     return prod
#
#
# print(factorialInt(200))
# print(factorialFloat(200))


# for the floating point version, 200! overflows, i.e. is larger than the largest possible float,
# but integer calculations are done to arbitrary precision

# Exercise 4.2
# def quadratic(a, b, c):
#     '''
#
#     :param a: quadratic coefficient
#     :param b: linear coeff
#     :param c: constant term
#     :return: prints out both solutions to the associated quadratic equation
#     '''
#     print('Solution one: ', (-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a))
#     print('Solution two: ', (-b - math.sqrt(b ** 2 - 4 * a * c)) / (2 * a))
#     print('Alternative form one:', (2 * c) / (-b - math.sqrt(b ** 2 - 4 * a * c)))
#     print('Alternative form one:', (2 * c) / (-b + math.sqrt(b ** 2 - 4 * a * c)))
#
#
# quadratic(0.001, 1000, 0.001)


# Exercise 4.3
# def f(x):
#     return x*(x - 1)
#
# delta = math.pow(10, -2)
# print( (f(1 + delta) - f(1) )/delta)
# for smaller and smaller delta, the approximation of the limit gets better, but at some point it gets worse due to
# numerical error arising from the difference of numbers that are nearly equal.
