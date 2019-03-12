from scipy import sqrt, sin, pi, linspace, array, arange, empty
from scipy.fftpack import rfft
from pylab import plot, show, xlabel, ylabel, imshow, gray

# Constants
grating_width = 200 * 10 ** -6  # m
W = 10 * grating_width
wavelength = 500 * 10 ** -9  # m
slit_width = 20 * 10 ** -6  # m
alpha = pi / slit_width
focal_length = 1  # nm
N = 10000  # number of sample points
screen_width = 0.1  # m
spacing = wavelength * focal_length / W  # m

def y(u):
    if abs(u) > grating_width / 2:
        return 0
    else:
        return sqrt(sin(alpha * u) ** 2)


sample_points_u = arange(0, N) * W / N - W / 2
y_points = array(list(map(y, sample_points_u)), float)
y_fourier = rfft(y_points)
intensity_profile = W ** 2 / N ** 2 * abs(y_fourier) ** 2

num_points = int(screen_width / spacing)
x_points = arange(-num_points, num_points) * spacing / 2

Iarray = empty([100, 2 * num_points], float)
for k in range(100):
    Iarray[k, 0 : num_points] = intensity_profile[0 : num_points][::-1]
    Iarray[k, num_points : 2 * num_points] = intensity_profile[0 : num_points]
imshow(Iarray)
gray()
show()


plot(x_points[num_points: 2 * num_points], intensity_profile[0: num_points], 'b')
plot(x_points[0: num_points], intensity_profile[0: num_points][::-1], 'b')
# plot(sample_points_u, y_points, 'o')
show()
