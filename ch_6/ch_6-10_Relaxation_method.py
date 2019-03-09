from scipy import exp

# Constants
accuracy = 10 ** -8

def f(c, x):
    return 1 - exp(- c * x)

def error(x1, x2, c):
    return (x1 - x2) / (1 - 1 / (c * exp(-c * x1)))

x1 = 0.5
x2 = f(2, x1)

while(abs(error(x1, x2, 2)) > accuracy):
    x1, x2 = x2, f(2, x2)

print(x2)
