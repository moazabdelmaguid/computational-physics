from scipy import exp, linspace
from pylab import plot, show, xlabel, ylabel
# Constants
accuracy = 10 ** -8

def g(c):
    iterations = 1
    def f(x):
        return 1 - exp(- c * x)

    def error(x1, x2):
        return (x1 - x2) / (1 - 1 / (c * exp(-c * x1)))

    x1 = 0.5 # starting value
    x2 = f(x1)
    while(abs(error(x1, x2)) > accuracy):
        x1, x2 = x2, f(x2)
        iterations += 1
    print('number of iterations = ', iterations)
    return x2

c = linspace(0, 3, 300)
x = list(map(g, c))
plot(c, x, 'o')
xlabel('c')
ylabel('x')
show()
