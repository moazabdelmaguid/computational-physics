from pylab import plot, show, ylim, xlabel, ylabel, imshow, gray, colorbar, scatter, hot
from numpy import linspace, sin, cos, loadtxt, pi, exp, sqrt, empty, arange, zeros, array, meshgrid, full, dot

#from vpython import sphere, vector

# Exercise 3.1
# sunspotData = loadtxt('../cpresources/sunspots.txt')
# x = sunspotData[0:1001,0]
# y = sunspotData[0:1001,1]
# plot(x, y)
# show()

# Exercise 3.2
# thetaValues = linspace(0, 2*pi, 50)
# x = 2*cos(thetaValues) + cos(2*thetaValues)
# y = 2*sin(thetaValues) - sin(2*thetaValues)
# plot(x, y)
# show()

# thetas = linspace(0, 10*pi, 500)
# r = exp(cos(thetas)) - 2*cos(4*thetas) + sin(thetas/12)**5
# x = r*cos(thetas)
# y = r*sin(thetas)
# plot(x, y)
# show()

# data = loadtxt("../cpresources/circular.txt",float)
# imshow(data, origin='lower')
# gray()
# colorbar()
# show()

# # Example 3.1
#
# wavelength = 5.0
# k = 2*pi/wavelength
# xi0 = 1.0
# separation = 20.0 # separation of centers in cm
# side = 100.0 # Side of the square in cm
# points = 500 # Number of grid points along each side
# spacing = side/points # Spacing of points in cm
#
# # Calculate the positions of the centers of the circles
# x1 = side/2 + separation/2
# y1 = side/2
# x2 = side/2 - separation/2
# y2 = side/2
#
# # Make an empty array to store the heights
# xi = empty([points, points], float)
#
# # Calculate the values in the array
# for i in range(points):
#     y = spacing*i
#     for j in range(points):
#         x = spacing*j
#         r1 = sqrt((x-x1)**2 + (y-y1)**2)
#         r2 = sqrt((x-x2)**2 + (y-y2)**2)
#         xi[i,j] = xi0*sin(k*r1) + xi0*sin(k*r2)
#
# # Make the plot
# imshow(xi,origin='lower',extent=[0,side,0,side])
# gray()
# show()

# Exercise 3.3
# heights = loadtxt('../cpresources/stm.txt',float)
# imshow(heights,origin='lower')
# gray()
# show()

# L = 5
# R = 0.3
# for i in range(-L,L+1):
#     for j in range(-L, L+1):
#         for k in range(-L,L+1):
#             sphere(pos=vector(i,j,k),radius = R)


# # Exercise 3.6
# rValues = arange(1, 4, 0.001)
# xArray = zeros(rValues.size, float) + 0.5
#
# #iterate 1000 times to settle down to a fixed point or limit cycle if it's going to
# for i in range(0,1000):
#     xArray = rValues*xArray*(1-xArray)
#
# #iterate another 1001 times to get one set of branches
# for i in range(0, 1001):
#     xArray = rValues * xArray * (1 - xArray)
#
# #plot these
# plot(rValues, xArray, '.', color="k")
# #iterate another 1001 times to get the even branches
# for i in range(0, 1001):
#     xArray = rValues * xArray * (1 - xArray)
# #plot the rest
# plot(rValues, xArray, '.', color="k")
# show()


# # Exercise 3.7 Mandelbrot set
N = 450 # number of grid points in each direction
numIterations = 150
grid = full([ N, N ], numIterations, int) # grid to plot

cx, cy = meshgrid(linspace(-2, 2, N), linspace(-2,2,N)) # Create grid of c values

z = zeros([N,N], complex)

for n in range(0, numIterations):
    z = z**2 + cx + cy*1j
    for i in range(0, N):
        for j in range(0, N):
            if abs(z[i,j]) > 2:
                grid[i,j] = n+1

imshow(grid,origin='lower',extent=[-2,2,-2,2])
hot()
colorbar()
show()

# Exercise 3.8
# data = loadtxt('../cpresources/millikan.txt', float)
# x = data[:,0]
# y = data[:,1]
# N = x.size
#
# Ex = 1/N*sum(x)
# Ey = 1/N*sum(y)
# Exx = 1/N*sum(x**2)
# Exy = 1/N*dot(x, y)
# m = (Exy - Ex*Ey)/(Exx - Ex**2)
# c = (Exx*Ey - Ex*Exy)/(Exx - Ex**2)
#
# newY = empty(N, float)
# for i in range(0, N):
#     newY[i] = m*x[i] + c
#
# plot(x, y, 'ko')
# plot(x, newY, 'k')
# show()















