from scipy import floor, linspace, array, zeros, copy
from scipy.fftpack import rfft, irfft
from pylab import plot, show, xlabel, ylabel

# square wave with amplitude 1 and frequency 1 Hz
def square_wave(t):
    if floor(2 * t) % 2 == 0:
        return 1
    else:
        return -1

# number of sample points in one cycle
N = 1000

sample_points = linspace(0, 1, N)
square_wave_samples = array(list(map(square_wave, sample_points)), float)
square_wave_fourier = rfft(square_wave_samples)
first_10_coefficients = zeros(len(square_wave_fourier), float)
first_10_coefficients[0 : 9] = copy(square_wave_fourier[0 : 9])
compressed_square_wave = irfft(first_10_coefficients)

plot(sample_points, square_wave_samples, 'k')
plot(sample_points, compressed_square_wave, 'r')
xlabel('t')
ylabel('f(t)')
show()

# Discarding too many high frequency modes leaves us with low frequency wiggles and a poorer approximation
