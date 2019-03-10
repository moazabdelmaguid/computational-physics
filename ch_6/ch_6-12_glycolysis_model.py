from scipy import sqrt
# Find the stationary points of the system of equations x'(t) = -x + a*y + x^2 * y, y'(t) = b - a*y - x^2 y
# Note that the exact solutions are x = b, y = b / (a + b^2)

# Constants
a = 1
b = 2
accuracy = 10 ** -6

# We use a relaxation method to solve x = y(a + x^2), y = b / (a + x^2):

def f1(x, y):
    return y * (a + x ** 2)


def g1(x, y):
    return b / (a + x ** 2)


def f2(x, y):
    return sqrt(b / y - a)


def g2(x, y):
    return x / (a + x ** 2)


def stationary_points(f, g):
    """

    :param f: function for x = f(x,y)
    :param g: function for y = g(x,y)
    :return: a list [x, y] of the fixed points of the above system of equations
    """
    iterations = 1

    def relative_change(x1, x2):
        return (x1 - x2) / x1

    # starting values
    x1 = 0.5
    y1 = 0.25
    x2 = f(x1, y1)
    y2 = g(x1, y1)

    while abs(relative_change(x1, x2)) > accuracy and abs(relative_change(y1, y2)) > accuracy:
        if iterations > 1000000:
            print('number of iterations exceeded maximum')
            return 'error'

        x1, x2, y1, y2 = x2, f(x2, y2), y2, g(x2, y2)
        iterations += 1
    print('number of iterations = ', iterations)
    return [ x2, y2 ]

# if we try to use x = f1(x,y) and y = g1(x,y) to solve for the stationary points, the relaxation method fails
# but if we use f2 and g2 as defined above, which are just rearrangements of the original equations,
# the method nowconverges
print(stationary_points(f2, g2))
