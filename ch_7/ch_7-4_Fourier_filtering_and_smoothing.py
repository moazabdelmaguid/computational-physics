from scipy import loadtxt, zeros, copy, floor
from scipy.fftpack import rfft, irfft
from pylab import plot, show, xlabel, ylabel

dow = loadtxt("../../cpresources/dow.txt", float)
# plot(dow)
# show()

dow_fourier = rfft(dow)
N = len(dow_fourier)
first_10_percent = zeros(N, float)
first_10_percent[0 : int(N / 10)] = copy(dow_fourier[0 : int(N / 10)])
dow_first_10_p = irfft(first_10_percent)
first_2_percent = zeros(N, float)
first_2_percent[0 : int(N / 50)] = copy(dow_fourier[0 : int(N / 50)])
dow_first_2_p =irfft(first_2_percent)

plot(dow, 'k')
plot(dow_first_10_p, 'b')
plot(dow_first_2_p, 'g')
show()
