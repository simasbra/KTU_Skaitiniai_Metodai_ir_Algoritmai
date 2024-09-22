# Niutono metodas

import sympy as sp
import numpy as np
import math
import matplotlib.pyplot as plt


def funk(x): return np.sin(x**2)+0.2
def dfunk(x): return 2*x*np.cos(x**2)


x0 = -0.5
xxx = np.linspace(-5, 5, 1000)
fff = funk(xxx)
plt.plot(xxx, fff, 'b-')
plt.plot([-5, 5], [0, 0], 'k--')
plt.plot(x0, 0, 'k*')

for it in range(100):
    f = funk(x0)
    df = dfunk(x0)
    x1 = x0-f/df
    plt.plot([x0, x0, x1], [0, f, 0], 'k-', linewidth=0.5)
    tksl = max([np.abs(x0-x1), np.abs(funk(x1))])
    print(tksl)
    print('it=', it, ' x1= ', x1, '  funk= ', funk(x1), 'tksl=', tksl)
    if (tksl < 1e-8):
        break
    x0 = x1

fff = funk(xxx)

plt.plot(x1, 0, 'rp')

# Teiloro eilute ir zeros funkcija

TE, x = sp.symbols(('TE', 'x'))
f = sp.sin(x)
print(f)

N = 15
x0 = 1
TE = f.subs(x, x0)
print(TE)
df = f
for i in range(1, N+1):
    df = df.diff(x)
    TE = TE+(x-x0)**i*df.subs(x, x0)/math.factorial(i)  # print(TE)

xxx = np.linspace(-8, 8, 500)
fff = np.zeros(np.size(xxx))
TTT = np.zeros(np.size(xxx))
for i in range(len(xxx)):
    fff[i] = f.subs(x, xxx[i])
    TTT[i] = TE.subs(x, xxx[i])
plt.plot(xxx, fff, 'b-')
plt.plot(xxx, TTT, 'r-')
plt.grid()
plt.ylim((-2, 2))

print(TE)
a = sp.Poly(TE, x)
print(a)
kf = np.array(a.all_coeffs())
print(kf)
saknys = np.roots(kf)

for i in range(len(saknys)):
    print(saknys[i])

for i in range(len(saknys)):
    if saknys.imag[i] == 0:
        plt.plot(saknys[i], 0, 'k*')
