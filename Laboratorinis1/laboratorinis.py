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
koeficientaiNeigX = [val * -1 if i % 2 != 0 else val for i, val in enumerate(koeficientai)]

# Iverciai
# Nustatykite daugianario 𝑓(𝑥𝑥šaknų intervalą, taikydami „grubų“ ir „tikslesnį“ įverčius.
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
# Grafiškai pavaizduokite daugianarį tokiame intervale, kad matytųsi abu įverčiai.
plotIverciai = numpy.arange(galutinis[0], galutinis[1] + 0.00001, 0.00001)
plotFunkcija = funkcijaG(plotIverciai)

pyplot.title("Funkcija g(x)")
pyplot.xlabel("x")
pyplot.ylabel("y")
pyplot.plot(plotIverciai, plotFunkcija, "m")
pyplot.xlim([-5, 5])
pyplot.axhline(y=0, lw=1, color='k')
# pyplot.ylim([-5, 5])
pyplot.grid()
grubusIntervalas = pyplot.scatter([-R_GRUBUS, R_GRUBUS], [0, 0], marker="v", color="b")
tiklsIntervalas = pyplot.scatter([-R_TIKSL_NEIG, R_TIKSL_TEIG], [0, 0], marker="*", color="r")

pyplot.legend((grubusIntervalas, tiklsIntervalas), ('"Grubaus" intervalo galai', '"Tikslesnio" intervalo galai'), scatterpoints=1, fontsize=8)
pyplot.show()
