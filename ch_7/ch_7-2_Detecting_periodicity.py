from scipy import loadtxt, exp, zeros, pi
from pylab import plot, show, xlabel, ylabel

sunspot_data = loadtxt("../../cpresources/sunspots.txt", float)
time = sunspot_data[:, 0]
sunspots = sunspot_data[:, 1]

## plot sunspot data
plot(time, sunspots)
xlabel('Number of months since Jan 1749')
ylabel('Number of sunspots')
show()

def dft(y):
    N = len(y)
    c = zeros(N // 2 + 1, complex)
    for k in range(N // 2 + 1):
        for n in range(N):
            c[k] += y[n] * exp(-2j * pi * k * n / N)
    return c


## Fourier transform sunspot data
fourier_data = dft(sunspots)

def mag_squared(a):
    return abs(a) ** 2


plot(list(map(mag_squared, fourier_data))[1 : 100])
xlabel('k')
ylabel('|c_k|^2')
show()

## The peak is roughly near k = 21, 22 which is associated with a sine wave of period N / k = 3143 / 21 = 145 months
## which agrees with the sunspot plot above
