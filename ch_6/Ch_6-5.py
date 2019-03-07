from numpy.linalg import solve
from numpy import array
from cmath import polar, phase

R1 = R3 = R5 = 1000 # resistance in ohms
R2 = R4 = R6 = 2000 # in ohms
C1 = 10 ** -6 # capacitance in farads
C2 = 0.5 * 10 ** -6 # in farads
xp = 3 # in volts
omega = 1000 # in hertz

A = array([ [ 1 / R1 + 1 / R4 + 1j * omega * C1, - 1j * omega * C1, 0 ],
            [ -1j * C1, 1 / R2 + 1 / R5 + 1j * omega * C1 + 1j * omega * C2, - 1j * omega * C2 ],
            [ 0, -1j * omega *C2, 1 / R3 + 1/ R6 + 1j * omega * C2 ]], complex)
v = array([ xp / R1, xp / R2, xp/ R3 ], complex)

x = solve(A, v)
for i in range(len(x)):
    r, theta = polar(x[i])
    print(r, 'V', theta, 'rad')
