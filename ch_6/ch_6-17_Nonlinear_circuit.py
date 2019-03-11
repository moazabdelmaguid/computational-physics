from scipy import exp, array, dot

# Constants
V_plus = 5  # Volts
R1 = 1000  # in ohms
R2 = 4000
R3 = 3000
R4 = 2000
I0 = 3 * 10 ** -9  # in amps
V_T = 0.05  # volts
accuracy = 0.0001  # volts


# Using relaxation method
# def f(v1, v2):
#     return v2 + V_T * log(1 - v1 / (I0 * R2) + (V_plus - v1) / (I0 * R1))
#     # return 1 / (1 / R1 + 1 / R2) * (V_plus / R1 - I0 * (exp((v1 - v2) / V_T) - 1))
#
#
# def g(v1, v2):
#     return R4 / R3 * (V_plus - v2) - R4 / R1 * (v1 - V_plus) - R4 / R2 * v1
#
#  # starting values
# x1 = 3.44694
# y1 = 2.82956  #8071712136
# x2 = f(x1, y1)
# y2 = g(x1, y1)
# print(abs(y1-y2))
# print(abs(x1-x2))
# while abs(x1 - x2) > accuracy or abs(y1 - y2) > accuracy:
#     x1, x2, y1, y2 = x2, f(x2, y2), y2, g(x2, y2)
#     print('x2 = ', x2)
#     print('y2 = ', y2)
#
# print(x2, y2)

## Can't get relaxation method to converge

# Using Newton's method
def f(arr):
    return array([ (arr[0] - V_plus) / R1 + arr[0] / R2 + I0 * (exp((arr[0] - arr[1]) / V_T) - 1) ,
                   (V_plus - arr[1]) / R3 - arr[1] / R4 + I0 * (exp((arr[0] - arr[1]) / V_T) - 1)
    ],float)

def grad_f(arr):
    return array([ [ 1 / R1 + 1 / R2 + I0 / V_T * exp((arr[0] - arr[1]) / V_T), -I0 / V_T * exp((arr[0] - arr[1]) / V_T) ],
                   [ I0 / V_T * exp((arr[0] - arr[1]) / V_T), -1 / R3 - 1 / R4 - I0 / V_T * exp((arr[0] - arr[1]) / V_T) ] ]
    , float)


def inverse_matrix(arr):
    return 1 / (arr[0, 0] * arr[1, 1] - arr[0, 1] * arr[1, 0]) * array([ [ arr[1, 1], -arr[0, 1] ], [ -arr[1, 0], arr[0, 0] ] ], float)


x1 = array([3, 2.4], float)
delta_x = dot(inverse_matrix(grad_f(x1)),f(x1))
while abs(delta_x[0]) > accuracy or abs(delta_x[1]) > accuracy:
    delta_x = dot(inverse_matrix(grad_f(x1)), f(x1))
    x1 -= delta_x

print('V1 = ', x1[0])
print('V2 = ', x1[1])

