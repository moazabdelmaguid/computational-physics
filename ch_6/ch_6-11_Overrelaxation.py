from scipy import exp

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

def g_overrelaxation(c, omega):
    iterations = 1

    def f(x):
        return 1 - exp(- c * x)

    def f_prime(x):
        return c * exp(- c * x)

    def error(x1, x2):
        return (x1 - x2) / (1 - 1 / ((1 + omega) * f_prime(x1) - omega))

    x1 = 0.5  # starting value
    x2 = (1 + omega) * f(x1) - omega * x1
    while abs(error(x1, x2)) > accuracy:
        x1, x2 = x2, (1 + omega) * f(x2) - omega * x2
        iterations += 1
    print('For omega = ', omega,'the number of iterations = ', iterations)
    return x2

print(g(2))
print(g_overrelaxation(2, 0.75))

# We might want to use an omega < 0 if at each iteration the function overshoots the solution, i.e., we start at x=2
# and the next iteration gives x=3, but the solution is 2.5
