from scipy.fftpack import fft2, ifft2
from scipy import loadtxt, exp, empty, real
from pylab import imshow, plot, show, gray
from numpy.fft import rfft2, irfft2

# Constants
sigma = 25
blurred_photo = loadtxt('../../cpresources/blur.txt', float)
y_dim, x_dim = blurred_photo.shape

def point_spread(x, y):
    return exp(- ( x ** 2 + y ** 2 ) / (2 * sigma ** 2))


# calculate point spread function for each point
point_spread_array = empty([ y_dim, x_dim ], float)
for i in range(y_dim):
    for j in range(x_dim):
        point_spread_array[i, j] = point_spread( (j + y_dim / 2) % y_dim - y_dim / 2, \
                                                 (i + x_dim / 2) % x_dim - x_dim / 2)

# Fourier transform both
blurred_photo_fourier = rfft2(blurred_photo)
point_spread_fourier = rfft2(point_spread_array)

# divide
unblurred_fourier = empty([ y_dim, x_dim // 2 + 1], complex)
epsilon = 10 ** -4
for i in range(x_dim // 2 + 1):
    for j in range(y_dim):
        if abs(point_spread_fourier[j, i]) < epsilon:
            unblurred_fourier[j, i] = blurred_photo_fourier[j, i]
        else:
            # Note we dont need a factor of x_dim * y_dim in the denominator due to normalization of rfft2
            unblurred_fourier[j, i] = blurred_photo_fourier[j, i] / (point_spread_fourier[j, i])

imshow(irfft2(unblurred_fourier))
gray()
show()

