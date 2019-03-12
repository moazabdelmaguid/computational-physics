from scipy.fftpack import fft
from scipy import loadtxt
from pylab import plot, show, xlabel, ylabel

piano = loadtxt('../../cpresources/piano.txt', float)
# plot(piano)
# show()

piano_fourier = fft(piano)
plot(abs(piano_fourier[0:9999]))
xlabel('k')
ylabel('|c_k|')
show()

trumpet = loadtxt('../../cpresources/trumpet.txt', float)
trumpet_fourier = fft(trumpet)
plot(abs(trumpet_fourier[0:9999]))
xlabel('k')
ylabel('|c_k|')
show()

# the frequency associated with k is given by f =  k * (44100 samples/s) / 100000 samples
