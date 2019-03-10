from scipy import tan, sqrt, linspace
from pylab import plot, show, xlabel

# Constants
w = 1  # width of well in nm
V = 20  # well depth in eV
m = 9.1094 * 10 ** -31  # electron mass in kg
h_bar = 6.582119514 * 10 ** -16  # h bar in eV * s
c = w ** 2 * m / (1.602 * 10 ** -1) / (2 * h_bar ** 2)
accuracy = 0.001  # in eV


def f(e):
    # use 1 eV = 1.602 * 10 ** -19 J to convert w^2 m in to ev * s^2
    return tan(sqrt(c * e))


def g(e):
    return sqrt((V - e) / e)


def h(e):
    return -sqrt(e / (V - e))


# make plots
E = linspace(0.05, 19.95, 100)
y1 = list(map(f, E))
y2 = list(map(g, E))
y3 = list(map(h, E))
plot(E, y1, 'o')
plot(E, y2, 'ko')
plot(E, y3, 'ro')
xlabel('E (eV)')
show()


def f1(x):
    return f(x) - g(x)


def f2(x):
    return f(x) - h(x)


# from exercise 6.13
def find_root(f, x1, x2, accuracy):
    """
    finds the root of f(x) = 0 in the interval (x1, x2)
    :param f: function f(x)
    :param x1: lower bound of interval
    :param x2: upper bound of interval
    :param accuracy: target accuracy
    :return: float
    """

    def midpoint(x, y):
        return (x + y) / 2

    def have_same_sign(x, y):
        if x < 0 and y < 0 or x > 0 and y > 0:
            return True
        else:
            return False

    #
    while abs(x1 - x2) > accuracy:
        x = midpoint(x1, x2)
        if have_same_sign(f(x1), f(x)):
            x1 = x
        elif have_same_sign(f(x), f(x2)):
            x2 = x
        elif abs(x) < accuracy:
            return x

    return midpoint(x1, x2)


print('E0 = ', find_root(f1, 0.02, 0.75, accuracy))
print('E1 = ', find_root(f2, 1, 1.5, accuracy))
print('E2 = ', find_root(f1, 2.5, 4.5, accuracy))
print('E3 = ', find_root(f2, 5, 6, accuracy))
print('E4 = ', find_root(f1, 7.5, 9, accuracy))
print('E5 = ', find_root(f2, 10, 12, accuracy))
