from scipy import floor, linspace, array, zeros, copy, loadtxt
from scipy.fftpack import rfft, irfft, dct, idct
from pylab import plot, show, xlabel, ylabel

dow2 = loadtxt('../../cpresources/dow2.txt', float)
# plot(dow2)
# show()

# Using the discrete Fourier transform
dow2_fourier = rfft(dow2)
N = len(dow2_fourier)
first_2_percent = zeros(N, float)
first_2_percent[0 : int(N / 50)] = copy(dow2_fourier[0 : int(N / 50)])
smoothed_dow2 = irfft(first_2_percent)

# using the discrete cosine transform
dow2_cos = dct(dow2)
n = len(dow2_cos)
first_cos_2_percent = zeros(n, float)
first_cos_2_percent[0 : int(n / 50)] = copy(dow2_cos[0 : int(n / 50)])
smoothed_cos_dow2 = idct(first_cos_2_percent) / (2*n) # need factor of 1 / 2n for normalization

plot(dow2, 'k')
plot(smoothed_dow2, 'g')
plot(smoothed_cos_dow2, 'r')
show()
