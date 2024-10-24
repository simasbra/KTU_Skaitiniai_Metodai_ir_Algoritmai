import numpy
import matplotlib.pyplot as pyplot
from atspindysSuQR import atspindys_su_qr
from niutonoMetodas import netiesiniu_lygciu_sistema

# 1 Tiesinių lygčių sistemų sprendimas

# Lentelėje 1 duotos tiesinės lygčių sistemos, 2 lentelėje nurodytas metodas ir lygčių sistemų
# numeriai (iš 1 lentelės). Reikia suprogramuoti nurodytą metodą ir išspręsti pateiktas lygčių
# sistemas. Programoje sprendžiant lygčių sistemas turi būti įvertinti atvejai:
# • kai lygčių sistema turi vieną sprendinį;
# • kai lygčių sistema sprendinių neturi;
# • kai lygčių sistema turi be gali daug sprendinių.
# Patikrinkite gautus sprendinius įrašydami juos į pradinę lygčių sistemą.
# Gautą sprendinį patikrinkite naudodami išorinius išteklius.

# Užduoties Nr.: 3
# Metodas: QR
# Lygčių sistemų Nr.: 7, 15, 16

A7 = numpy.matrix([[5, 1, 3, -1],
                   [1, 4, 0, 1],
                   [3, 0, 11, 5],
                   [-1, 2, 5, 2]]).astype(float)
b7 = numpy.matrix([-5, 3, -4, -7]).transpose().astype(float)

A15 = numpy.matrix([[-4, 3, -5, 5],
                    [2, 2, 4, 3],
                    [0, -7, 2, -5],
                    [1, 1, 2, 1]]).astype(float)
b15 = numpy.matrix([-4, 16, 9, 8]).transpose().astype(float)

A16 = numpy.matrix([[1, 2, 1, 1],
                    [2, -5, 1, 2],
                    [4, -1, 3, 4],
                    [3, -3, 2, 3]]).astype(float)
b16 = numpy.matrix([-7, 3, 0, 2]).transpose().astype(float)

print("\nA7\n")
atspindys_su_qr(A7, b7)

print("\nA15\n")
atspindys_su_qr(A15, b15)

print("\nA16\n")
atspindys_su_qr(A16, b16)

# 2 Netiesinių lygčių sistemų sprendimas

# Duota netiesinių lygčių sistema (3 lentelė):
# 𝑍1(𝑥1, 𝑥2) = 0
# 𝑍2(𝑥1, 𝑥2) = 0
# a. Skirtinguose grafikuose pavaizduokite paviršius 𝑍1(𝑥1, 𝑥2) ir 𝑍2(𝑥1, 𝑥2).
# b. Užduotyje pateiktą netiesinių lygčių sistemą išspręskite grafiniu būdu.
# c. Nagrinėjamoje srityje sudarykite stačiakampį tinklelį (𝑥1, 𝑥2 poras). Naudodami užduotyje nurodytą metodą
# apskaičiuokite netiesinių lygčių sistemos sprendinius, kai pradinis artinys įgyja tinklelio koordinačių reikšmes.
# Tinklelyje vienodai pažymėkite taškus, kuriuos naudojant kaip pradinius artinius gaunamas tas pats sprendinys.
# Lentelėje pateikite apskaičiuotus skirtingus sistemos prendinius ir bent po vieną jam atitinkantį pradinį artinį.
# d. Gautus sprendinius patikrinkite naudodami išorinius išteklius

# Lygčių sistema:
# 𝑥1^2 + (𝑥2 + cos(𝑥1))^2 − 40 = 0
# (𝑥1 / 2 )^3 + 25 * 𝑥2^2 − 50 = 0
# Metodas: Niutono

# 3a. Skirtinguose grafikuose pavaizduokite paviršius 𝑍1(𝑥1, 𝑥2) ir 𝑍2(𝑥1, 𝑥2).

xNTLS1 = numpy.linspace(-5, 5, 200)
xNTLS2 = numpy.linspace(-5, 5, 200)
X1, X2 = numpy.meshgrid(xNTLS1, xNTLS2)

Z1 = X1**2 + (X2 + numpy.cos(X1))**2 - 40
Z2 = (X1 / 2)**3 + 25 * X2**2 - 50

# Z1
figureZ1 = pyplot.figure()
axeZ1 = pyplot.axes(projection="3d")
surfaceZ1 = axeZ1.plot_surface(X1, X2, Z1, cmap=pyplot.cm.plasma, alpha=0.4)
contourZ1 = axeZ1.contour(X1, X2, Z1, levels=[0], colors="r", linestyles="solid", offset=0)
axeZ1.set_xlabel("x", labelpad=20)
axeZ1.set_ylabel("y", labelpad=20)
axeZ1.set_zlabel("z", labelpad=20)

pyplot.show()

# Z2
figureZ2 = pyplot.figure()
axeZ2 = pyplot.axes(projection="3d")
surfaceZ2 = axeZ2.plot_surface(X1, X2, Z2, cmap=pyplot.cm.plasma, alpha=0.4)
contourZ2 = axeZ2.contour(X1, X2, Z2, levels=[0], colors="r", linestyles="solid", offset=0)
axeZ2.set_xlabel("x", labelpad=20)
axeZ2.set_ylabel("y", labelpad=20)
axeZ2.set_zlabel("z", labelpad=20)

pyplot.show()

# 3b. Užduotyje pateiktą netiesinių lygčių sistemą išspręskite grafiniu būdu.

figureNTLS = pyplot.figure(1, figsize=pyplot.figaspect(0.45))
axeNTLS1 = figureNTLS.add_subplot(1, 2, 1, projection="3d")
axeNTLS1.set_xlabel("x")
axeNTLS1.set_ylabel("y")
axeNTLS1.set_zlabel("z")

axeNTLS2 = figureNTLS.add_subplot(1, 2, 2)
axeNTLS2.set_xlabel("x")
axeNTLS2.set_ylabel("y")
axeNTLS2.set_title("Grafinis netiesinių TLS sprendimas")

pyplot.draw()

xNTLS1 = numpy.linspace(-8, 8, 50)
xNTLS2 = numpy.linspace(-8, 8, 50)
nuliaiNTLS = numpy.zeros(shape=(len(xNTLS1), len(xNTLS2), 2))

X1, X2 = numpy.meshgrid(xNTLS1, xNTLS2)
for i in range(0, len(xNTLS1)):
    for j in range(0, len(xNTLS2)):
        nuliaiNTLS[i, j, :] = netiesiniu_lygciu_sistema([X1[i][j], X2[i][j]]).transpose()

surf1 = axeNTLS1.plot_surface(X1, X2, nuliaiNTLS[:, :, 0], color="blue", alpha=0.4, linewidth=0.1, antialiased=True)
CS11 = axeNTLS1.contour(X1, X2, nuliaiNTLS[:, :, 0], [0], colors="b")

surf2 = axeNTLS1.plot_surface(X1, X2, nuliaiNTLS[:, :, 1], color="red", alpha=0.4, linewidth=0.1, antialiased=True)
CS12 = axeNTLS1.contour(X1, X2, nuliaiNTLS[:, :, 1], [0], colors="r")

CS1 = axeNTLS2.contour(X1, X2, nuliaiNTLS[:, :, 0], [0], colors="b")
CS2 = axeNTLS2.contour(X1, X2, nuliaiNTLS[:, :, 1], [0], colors="r")

pyplot.grid()
pyplot.show()

# c. Nagrinėjamoje srityje sudarykite stačiakampį tinklelį (𝑥1, 𝑥2 poras). Naudodami užduotyje nurodytą metodą
# apskaičiuokite netiesinių lygčių sistemos sprendinius, kai pradinis artinys įgyja tinklelio koordinačių reikšmes.
# Tinklelyje vienodai pažymėkite taškus, kuriuos naudojant kaip pradinius artinius gaunamas tas pats sprendinys.
# Lentelėje pateikite apskaičiuotus skirtingus sistemos prendinius ir bent po vieną jam atitinkantį pradinį artinį.

# ......
