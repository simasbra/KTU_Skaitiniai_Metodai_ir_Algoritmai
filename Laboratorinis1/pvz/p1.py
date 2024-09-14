# funkcijos grafiko vaizdavimas
import sympy as sp
import numpy as np
import math
import matplotlib.pyplot as plt


def funk(x): return np.sin(x)


Pi = math.pi
xxx = np.linspace(-2*Pi, 2*Pi, 50)  # print(xxx)
fff = funk(xxx)
print(fff)
plt.plot(xxx, fff, 'b-')
plt.grid()
plt.xlabel('x')
plt.ylabel('f')

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
f1, = ax1.plot(xxx, fff, 'b-.')
ax1.grid()
ax1.set_xlabel('x')
ax1.set_ylabel('f')
f2, = ax2.plot(xxx, fff-xxx, 'r-*')
ax2.grid()
ax2.set_xlabel('x')
ax2.set_ylabel('f')


# simboline israiska duotos funkcijos grafiko vaizdavimas


x, f = sp.symbols(('x', 'f'))
print(x, f)
f = sp.sin((x**2))
print(f)
Pi = math.pi
xxx = np.linspace(-2*Pi, 2*Pi, 500)  # print(xxx)
fff = np.array(np.zeros(np.size(xxx)))  # print(fff)
for i in range(len(xxx)):
    fff[i] = f.subs(x, xxx[i]).evalf()

plt.plot(xxx, fff, 'b-')
plt.grid()
plt.xlabel('x')
plt.ylabel('f')
plt.show()

df = sp.diff(f, x)  # print(df)
dfff = np.array(np.zeros(np.size(xxx)))  # print(dfff)
for i in range(len(xxx)):
    dfff[i] = df.subs(x, xxx[i]).evalf()
# print(dfff)
plt.plot(xxx, dfff, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('df')


# saknu atskyrimo metodai: pusiaukirta


def funk(x): return np.sin(x)


Pi = math.pi
xn = -6
xn1 = 2

xxx = np.linspace(-2*Pi, 2*Pi, 50)  # print(xxx)
fff = funk(xxx)

tikslumas = 9999
it = 0
while tikslumas > 1e-6:
    it += 1
    fn = funk(xn)
    fn1 = funk(xn1)
    xmid = (xn+xn1)/2
    fmid = funk(xmid)  # pusiaukirta

    plt.plot(xxx, fff, 'b-.')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('f')
    plt.plot(xn, 0, 'k*')
    plt.plot(xn1, 0, 'b*')
    plt.plot(xmid, 0, 'rp')
    plt.plot([xn, xn], [0, fn], 'k--')
    plt.plot([xn1, xn1], [0, fn1], 'k--')

    if np.sign(fn) == np.sign(fmid):
        xn = xmid
    else:
        xn1 = xmid
    tikslumas = np.abs(xn-xn1)
    print('iteracija ', it, 'saknis= ', xmid, ' tikslumas= ', tikslumas)
    plt.show()
