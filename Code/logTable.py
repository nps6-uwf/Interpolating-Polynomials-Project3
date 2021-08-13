from math import log10
from scipy.interpolate import lagrange

x = [2, 3]
y = [log10(i) for i in x]
poly = lagrange(x, y) # Lagrange model

print(poly(2.5), log10(2.5), abs(poly(2.5) - log10(2.5)))
