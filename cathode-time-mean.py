from math import *

#chisq = 6.52543
#p0 = 0.0295522 +- 0.000221807
#p1 = -0.356829 +- 0.167502
#p2 = 3.04875 +- 0.18189
#p3 = 51.3332 +- 0.208807
#p4 = 58.5917 +- 0.0439822

A = 0.0367744
a = -0.0501774
b = 2.74145
c = 51.3882
d = 59.0447

# integral = int (f(x) * x) dx / int (f(x)) dx

i0 = 0.5 * A * (b - a)
i1 = A * (c - b)
i2 = 0.5 * A * (d - c)
print(i0, i1, i2)
m0 = (A / (b - a)) * ((1.0 / 3.0) * (b**3.0 - a**3.0) + 0.5 * a**3.0 - 0.5 * a * b**2.0)
m1 = 0.5 * A * (c**2.0 - b**2.0)
m2 = 0.5 * A * (d**2.0 - c**2.0) - (A / (d - c)) * ((1.0/3.0) * (d**3.0 - c**3.0) - 0.5 * (c * d**2.0 - c**3.0))
print(m0, m1, m2)
integral = i0 + i1 + i2
mean = m0 + m1 + m2
print(integral)
print(mean)
mean /= integral
print(mean)
mean1 = mean

A = 0.0607391
a = 2.08141
b = 5.50607
c = 35.708
d = 37.2706

i0 = 0.5 * A * (b - a)
i1 = A * (c - b)
i2 = 0.5 * A * (d - c)
print(i0, i1, i2)
m0 = (A / (b - a)) * ((1.0 / 3.0) * (b**3.0 - a**3.0) + 0.5 * a**3.0 - 0.5 * a * b**2.0)
m1 = 0.5 * A * (c**2.0 - b**2.0)
m2 = 0.5 * A * (d**2.0 - c**2.0) - (A / (d - c)) * ((1.0/3.0) * (d**3.0 - c**3.0) - 0.5 * (c * d**2.0 - c**3.0))
print(m0, m1, m2)
integral = i0 + i1 + i2
mean = m0 + m1 + m2
print(integral)
print(mean)
mean /= integral
print(mean)
mean2 = mean

print(mean2 / mean1)