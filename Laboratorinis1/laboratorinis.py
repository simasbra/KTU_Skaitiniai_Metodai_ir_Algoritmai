import sympy
import numpy
import math
import matplotlib.pyplot as pyplot


def daugianarisF(x):
    return -0.45 * x**4 + 1.04 * x**3 + 1.42 * x**2 - 2.67 * x - 0.97


def funkcijaG(x):
    return (numpy.cos(2 * x) / (numpy.sin(x) + 1.5)) - (x / 5)


def grubusIntervalas(koeficientai):
    max_koef = max(koeficientai[:-1], key=abs)
    R = 1 + (max_koef / koeficientai[-1])
    return round(R, 2)


def tikslesnisIntervalas(koeficientai):
    n = len(koeficientai)
    max_abs_index, B = min(enumerate(koeficientai[:-1]), key=lambda x: x[1])
    k = n - max_abs_index

    R = 1 + (abs(B) / koeficientai[-1])**(1 / k)
    return round(R, 2)


koeficientai = [-0.45, 1.04, 1.42, -2.67, -0.97]
koeficientaiNeigX = [
    val * -1 if i % 2 != 0 else val for i, val in enumerate(koeficientai)]

# 1.1
# Iverciai
# Nustatykite daugianario 𝑓(𝑥)šaknų intervalą,
# taikydami „grubų“ ir „tikslesnį“ įverčius.
if koeficientai[-1] < 0:
    koeficientai = [x * (-1) for x in koeficientai]

if koeficientaiNeigX[-1] < 0:
    koeficientaiNeigX = [x * (-1) for x in koeficientaiNeigX]

R_GRUBUS = grubusIntervalas(koeficientai)
grubus = (-R_GRUBUS, R_GRUBUS)
print(f'"grubus" intervalas: {grubus}')

R_TIKSL_TEIG = tikslesnisIntervalas(koeficientai)
R_TIKSL_NEIG = tikslesnisIntervalas(koeficientaiNeigX)
tikslesnis = (-R_TIKSL_NEIG, R_TIKSL_TEIG)
print(f'"tikslesnis" intervalas: {tikslesnis}')

galutinis = (-min(R_GRUBUS, R_TIKSL_NEIG), min(R_GRUBUS, R_TIKSL_TEIG))
print(f'galutinis intervalas: {galutinis}')

# Grafinis pavaizdavimas
# Grafiškai pavaizduokite daugianarį tokiame intervale,
# kad matytųsi abu įverčiai.
dfx = 0.00001
iverciaiF = numpy.arange(galutinis[0], galutinis[1] + dfx, dfx)
figureF = daugianarisF(iverciaiF)
plotF, (plotF1, plotF2) = pyplot.subplots(1, 2, figsize=(10, 5))

# pirmojo grafiko braizymas
plotF1.set_title("-0.45x^4 + 1.04x^3 + 1.42x^2 - 2.67x - 0.97")
plotF1.set_xlabel("x")
plotF1.set_ylabel("y")
plotF1.plot(iverciaiF, figureF, "m")
plotF1.axhline(y=0, lw=1, color='k')
plotF1.grid()
grubusIntervalas = plotF1.scatter(
    [-R_GRUBUS, R_GRUBUS], [0, 0], marker="v", color="b")
tiklsIntervalas = plotF1.scatter(
    [-R_TIKSL_NEIG, R_TIKSL_TEIG], [0, 0], marker="*", color="r")

plotF1.legend((grubusIntervalas, tiklsIntervalas),
              ('"Grubaus" intervalo galai', '"Tikslesnio" intervalo galai'),
              scatterpoints=1, fontsize=8)

# antro grafiko braizymas aiskiam vizualizavimui
plotF2.set_title(
    "-0.45x^4 + 1.04x^3 + 1.42x^2 - 2.67x - 0.97\npriartintas vaizdas geresniam vizualizavimui")
plotF2.set_xlabel("x")
plotF2.set_ylabel("y")
plotF2.plot(iverciaiF, figureF, "m")
plotF2.set_xlim([-4, 4])
plotF2.set_ylim([-2, 2])
plotF2.axhline(y=0, lw=1, color='k')
plotF2.grid()
grubusIntervalas = plotF2.scatter(
    [-R_GRUBUS, R_GRUBUS], [0, 0], marker="v", color="b")
tiklsIntervalas = plotF2.scatter(
    [-R_TIKSL_NEIG, R_TIKSL_TEIG], [0, 0], marker="*", color="r")

plotF2.legend((grubusIntervalas, tiklsIntervalas),
              ('"Grubaus" intervalo galai', '"Tikslesnio" intervalo galai'),
              scatterpoints=1, fontsize=8)
pyplot.show()

# Funkciją 𝑔(𝑥) grafišksi pavaizduokite užduotyje nurodytame intervale.
dgx = 0.00001
intervalaiG = numpy.arange(-5, 5+dgx, dgx)
figureG = funkcijaG(intervalaiG)
plotG, (plotG1, plotG2) = pyplot.subplots(1, 2, figsize=(10, 5))

# pirmojo grafiko braizymas pagal salygas
plotG1.set_title(
    "cos(2x) / (sin(x) + 1.5)) - (x / 5)\nkai -5 <= x <= 5")
plotG1.set_xlabel("x")
plotG1.set_ylabel("y")
plotG1.plot(intervalaiG, figureG, 'r')
plotG1.axhline(y=0, lw=1, color='k')
plotG1.grid()

# antro grafiko braizymas aiskiam saknu vizualizavimui
plotG2.set_title(
    "cos(2x) / (sin(x) + 1.5)) - (x / 5)\npriartintas vaizdas geresniam vizualizavimui")
plotG2.set_xlabel("x")
plotG2.set_ylabel("y")
plotG2.plot(intervalaiG, figureG, 'b')
plotG2.axhline(y=0, lw=1, color='k')
plotG2.set_xlim([-3, 4])
plotG2.set_ylim([-2, 2])
plotG2.grid()

pyplot.tight_layout()
pyplot.show()
