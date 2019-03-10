from scipy import exp

# part a: solve 5 * e^-x + x - 5 = 0 using binary search
# Constants
accuracy = 10 ** -6

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

def f(x):
    return 5 * exp(-x) + x - 5

print(find_root(f, 4, 6, accuracy))
