from scipy import exp

# Constants
accuracy = 10 ** -8

def g(c):

    def f(x):
        return 1 - exp(- c * x)

    def error(x1, x2):
        return (x1 - x2) / (1 - 1 / (c * exp(-c * x1)))

    x1 = 0.5 # starting value
    x2 = f(x1)
    while(abs(error(x1, x2, c)) > accuracy):
        x1, x2 = x2, f(c, x2)

    return x2


