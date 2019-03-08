from vpython import sphere, vector, rate
from scipy import sin, cos, pi, arange, array







s = sphere(pos = vector( 1, 0, 0 ), radius = 0.1)
for theta in arange(0, 10 * pi, 0.1):
    rate(30)
    x = cos(theta)
    y = sin(theta)
    s.pos = vector(x, y, 0)
