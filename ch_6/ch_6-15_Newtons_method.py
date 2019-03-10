from scipy import linspace
from pylab import plot, show

def f(x):
    return 924 * x ** 6 - 2772 * x ** 5 + 3150 * x **  4 - 1680 * x ** 3 + 420 * x ** 2 - 42 * x + 1


def f_prime(x):
    return 5544 * x ** 5 - 13860 * x ** 4 + 12600 * x ** 3 - 5040 * x ** 2 + 840 * x - 42

# Plot f from x = 0 to 1
x = linspace(0, 1, 100)
f_vals = list(map(f, x))
plot(x, f_vals, 'o')
show()

# Solve for the roots of f using Newton's method
accuracy = 10 ** -10

def find_root(f, f_prime, start_val, accuracy):
    x = start_val
    delta = f(x) / f_prime(x)
    while abs(delta) > accuracy:
        delta = f(x) / g(x)
        x -= delta

    return x

print('The six roots of P_6(x) are')
print(find_root(f, f_prime, 0, accuracy))
print(find_root(f, f_prime, 0.2, accuracy))
print(find_root(f, f_prime, 0.4, accuracy))
print(find_root(f, f_prime, 0.6, accuracy))
print(find_root(f, f_prime, 0.8, accuracy))
print(find_root(f, f_prime, 0.95, accuracy))

