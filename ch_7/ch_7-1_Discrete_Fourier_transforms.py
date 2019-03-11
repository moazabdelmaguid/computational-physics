from scipy import exp, pi, zeros, empty, sin
from pylab import plot, show

# Constants
N = 1000  # number of sample points

def dft(y):
    N = len(y)
    c = zeros(N // 2 + 1, complex)
    for k in range(N // 2 + 1):
        for n in range(N):
            c[k] += y[n] * exp(-2j * pi * k * n / N)
    return c


# square wave with amplitude 1
# create the sample points
# square_wave = zeros(1000, complex)
# for i in range(1, N // 2):
#     square_wave[i] = 1.0
# # print(list(map(abs,dft(square_wave))))
# plot(list(map(abs,dft(square_wave))), 'o')
# show()

# sawtooth wave y_n = n
# sawtooth = empty(N, float)
# for i in range(N):
#     sawtooth[i] = i
# # plot(sawtooth)
# plot(list(map(abs,dft(sawtooth))), 'o')
# show()

# Modulated sine wave
mod_sine = empty(N, float)
for n in range(N):
    mod_sine[n] = sin(pi * n / N) * sin(20 * pi * n / N)

plot(list(map(abs,dft(mod_sine))), 'o')
show()
