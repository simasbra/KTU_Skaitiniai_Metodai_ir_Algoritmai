import sympy
import numpy
import math
import matplotlib.pyplot as pyplot
from scipy.optimize import fsolve


# -0.45x^4 + 1.04x^3 + 1.42x^2 - 2.67x - 0.97
def daugianarisF(x):
    return -0.45 * x**4 + 1.04 * x**3 + 1.42 * x**2 - 2.67 * x - 0.97


# cos(2x) / (sin(x) + 1.5)) - (x / 5), kai -5 <= x <= 5
def funkcijaG(x):
    return (numpy.cos(2 * x) / (numpy.sin(x) + 1.5)) - (x / 5)


def grubus_intervalas(koeficientai):
    max_koef = max(koeficientai[:-1], key=abs)
    R = 1 + (max_koef / koeficientai[-1])
    return round(R, 2)


def tikslesnis_intervalas(koeficientai):
    n = len(koeficientai)
    max_abs_index, B = min(enumerate(koeficientai[:-1]), key=lambda x: x[1])
    k = n - max_abs_index
    R = 1 + (abs(B) / koeficientai[-1])**(1 / k)
    return round(R, 2)


def skenavimas(funkcija, intervalas, zingsnis):
    # pradinė x reikšmė - intervalo pradžia
    x = intervalas[0]
    # čia bus laikomi šaknų atskyrimo intervalai
    saknys = []
    # nusprendžiama, kurią funkciją naudosime (pagal funkcijos parametrą)
    func = daugianarisF if funkcija == 'f' else funkcijaG
    # nustatomas funkcijos ženklas intervalo pradžioje
    dabZenklas = False if func(x) < 0 else True
    # kol prieiname intervalo galą tikriname ženklus tarp skirtingų galų
    # jeigu ženklai nesutampa - tame intervale bus šaknis, tad ją ir išsaugome
    iterations = 0
    while x < intervalas[1]:
        buvZenklas = dabZenklas
        dabZenklas = False if func(x) < 0 else True

        if buvZenklas != dabZenklas:
            saknys.append((x-zingsnis, x))
        x += zingsnis
        iterations += 1
    return saknys, iterations


def pusiaukirtos(funkcija, saknys, tolerancija=1e-10, iteracijosMax=1000):
    patikslintosSaknys = []
    iterations = []
    for saknuPora in saknys:
        xFrom = saknuPora[0]
        xTo = saknuPora[1]
        xMid = (xFrom + xTo) / 2
        iteration = 0
        while (numpy.abs(funkcija(xMid)) > tolerancija):
            iteration += 1
            if (numpy.sign(funkcija(xMid)) == numpy.sign(funkcija(xFrom))):
                xFrom = xMid
            else:
                xTo = xMid
            xMid = (xFrom + xTo) / 2
        patikslintosSaknys.append(xMid)
        iterations.append(iteration)
        # print(f'Artinys ({xFrom}, {xTo})')
        # print(f'Patikslinta saknis: {xMid} Iteracijos: {format(iteration)}')
    return patikslintosSaknys, iterations


def kvazi_niutono(funkcija, saknys, iteracijosMax=1000):
    patikslintosSaknys = []
    iterations = []
    for saknuPora in saknys:
        x0 = saknuPora[0]
        x1 = saknuPora[1]
        iteration = 0
        for i in range(iteracijosMax):
            funkcijax0 = funkcija(x0)
            funkcijax1 = funkcija(x1)
            if funkcijax1 == funkcijax0:  # Patikrina, ar yra dalyba iš nulio
                break
            x2 = x1 - (funkcija(x1) * (x1 - x0) /
                       float(funkcijax1 - funkcijax0))
            x0, x1 = x1, x2
            iteration += 1
        patikslintosSaknys.append(x2)
        iterations.append(iteration)
    return patikslintosSaknys, iterations


koeficientai = [-0.45, 1.04, 1.42, -2.67, -0.97]
koeficientaiNeigX = [
    val * -1 if i % 2 != 0 else val for i, val in enumerate(koeficientai)]

# 1.1
# Iverciai
# Nustatykite daugianario 𝑓(𝑥)šaknų intervalą,
# taikydami „grubų“ ir „tikslesnį“ įverčius.

print("\n1.1")
if koeficientai[-1] < 0:
    koeficientai = [x * (-1) for x in koeficientai]

if koeficientaiNeigX[-1] < 0:
    koeficientaiNeigX = [x * (-1) for x in koeficientaiNeigX]

R_GRUBUS = grubus_intervalas(koeficientai)
grubus = (-R_GRUBUS, R_GRUBUS)
print(f'"grubus" intervalas: {grubus}')

R_TIKSL_TEIG = tikslesnis_intervalas(koeficientai)
R_TIKSL_NEIG = tikslesnis_intervalas(koeficientaiNeigX)
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
grubus_intervalas = plotF1.scatter(
    [-R_GRUBUS, R_GRUBUS], [0, 0], marker="v", color="b")
tiklsIntervalas = plotF1.scatter(
    [-R_TIKSL_NEIG, R_TIKSL_TEIG], [0, 0], marker="*", color="r")

plotF1.legend((grubus_intervalas, tiklsIntervalas),
              ('"Grubaus" intervalo galai', '"Tikslesnio" intervalo galai'),
              scatterpoints=1, fontsize=8)

# antro grafiko braizymas aiskiam vizualizavimui
plotF2.set_title(
    "-0.45x^4 + 1.04x^3 + 1.42x^2 - 2.67x - 0.97\npriartintas vaizdas geresniam vizualizavimui")
plotF2.set_xlabel("x")
plotF2.set_ylabel("y")
plotF2.plot(iverciaiF, figureF, "m")
plotF2.set_xlim([-2, 2])
plotF2.set_ylim([-2, 2])
plotF2.axhline(y=0, lw=1, color='k')
plotF2.grid()
grubus_intervalas = plotF2.scatter(
    [-R_GRUBUS, R_GRUBUS], [0, 0], marker="v", color="b")
tiklsIntervalas = plotF2.scatter(
    [-R_TIKSL_NEIG, R_TIKSL_TEIG], [0, 0], marker="*", color="r")

plotF2.legend((grubus_intervalas, tiklsIntervalas),
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

# 1.2
# Naudodami skenavimo algoritmą su nekintančiu skenavimo žingsniu
# raskite šaknų atskyrimo intervalus.

print("\n1.2")
zingsnis = 0.1
daugianarioSaknys, iterationsSkenF = skenavimas('f', galutinis, zingsnis)
funkcijosSaknys, iterationsSkenG = skenavimas('g', (-5, 5), zingsnis)
print(f'Daugianario saknys skenavimo algoritmu: {daugianarioSaknys}')
print(f'Funkcijos saknys skenavimo algoritmu: {funkcijosSaknys}')

# 1.3
# Skenavimo metodu atskirtas daugianario ir funkcijos šaknis tikslinkite
# užduotyje nurodytais metodais. Skaičiavimo scenarijuje turi būti
# panaudotos skaičiavimų pabaigos sąlygos.

print("\n1.3")
# pusiaukirtos metodas
epsilonF = 1e-15
pusiaukirtosF, iterationsPusF = pusiaukirtos(
    daugianarisF, daugianarioSaknys, epsilonF, 1000)
# print(f'Daugianario saknys pusiaukirtos metodu: {pusiaukirtosF}')

epsilonG = 1e-14
pusiaukirtosG, iterationsPusG = pusiaukirtos(
    funkcijaG, funkcijosSaknys, epsilonG, 1000)
# print(f'Funkcijos saknys pusiaukirtos metodu: {pusiaukirtosG}')

# kvazi-niutono (kirstiniu) metodas
kvaziNiutonoF, iterationsKvazF = kvazi_niutono(
    daugianarisF, daugianarioSaknys, 1000)
# print(f'Daugianario saknys Kvazi-Niutono metodu: {kvaziNiutonoF}')

kvaziNiutonoG, iterationsKvazG = kvazi_niutono(
    funkcijaG, funkcijosSaknys, 1000)
# print(f'Funkcijos saknys Kvazi-Niutono metodu: {kvaziNiutonoG}')

# Skaičiavimų rezultatus pateikite lentelėje, kurioje nurodykite šaknies
# tikslinimui naudojamą metodą, pradinį artinį arba atskyrimo intervalą, gautą
# sprendinį (šaknį), funkcijos reikšmę ties šaknimi, tikslumą, iteracijų skaičių.
# Palyginkite, kuriuo metodu sprendiniui rasti panaudota mažiau iteracijų;
print('-' * 107)
print("| METODAS".ljust(15) + " | ATSKYRIMO INTERVALAS".ljust(25) +
      " | GAUTA ŠAKNIS".ljust(20) + " | F. REIKŠMĖ".ljust(20) +
      " | epsilon".ljust(10) + " | ITERACIJOS".ljust(15) + " |")
print('-' * 107)
print("| f(x) = -0.45x^4 + 1.04x^3 + 1.42x^2 - 2.67x - 0.97" + ' ' * 54 + '|')
print('-' * 107)
for i in range(len(daugianarioSaknys)):
    saknuPora = daugianarioSaknys[i]
    pusiaukirtos = pusiaukirtosF[i]
    iterationPus = iterationsPusF[i]
    kvaziNiuntono = kvaziNiutonoF[i]
    iterationKvaz = iterationsKvazF[i]
    print("| Pusiaukirtos".ljust(15) +
          f" | [{round(saknuPora[0], 3)} - {round(saknuPora[1], 3)}]".ljust(25) +
          f" | {round(pusiaukirtos, 10)}".ljust(20) +
          f" | {round(daugianarisF(pusiaukirtos), 14)}".ljust(20) +
          f" | {epsilonF}".ljust(10) +
          f" | {iterationPus}".ljust(15) + " |")
    print("| Kvazi-Niutono".ljust(15) +
          f" | [{round(saknuPora[0], 3)} - {round(saknuPora[1], 3)}]".ljust(25) +
          f" | {round(kvaziNiuntono, 10)}".ljust(20) +
          f" | {round(daugianarisF(kvaziNiuntono), 14)}".ljust(20) +
          f" | {epsilonF}".ljust(10) +
          f" | {iterationKvaz}".ljust(15) + " |")
print('-' * 107)
print("| g(x) = cos(2x) / (sin(x) + 1.5)) - (x / 5), kai -5 <= x <= 5" +
      ' ' * 44 + '|')
print('-' * 107)
for i in range(len(funkcijosSaknys)):
    saknuPora = funkcijosSaknys[i]
    pusiaukirtos = pusiaukirtosG[i]
    iterationPus = iterationsPusG[i]
    kvaziNiuntono = kvaziNiutonoG[i]
    iterationKvaz = iterationsKvazG[i]
    print("| Pusiaukirtos".ljust(15) +
          f" | [{round(saknuPora[0], 3)} - {round(saknuPora[1], 3)}]".ljust(25) +
          f" | {round(pusiaukirtos, 10)}".ljust(20) +
          f" | {round(daugianarisF(pusiaukirtos), 13)}".ljust(20) +
          f" | {epsilonG}".ljust(10) +
          f" | {iterationPus}".ljust(15) + " |")
    print("| Kvazi-Niutono".ljust(15) +
          f" | [{round(saknuPora[0], 3)} - {round(saknuPora[1], 3)}]".ljust(25) +
          f" | {round(kvaziNiuntono, 10)}".ljust(20) +
          f" | {round(daugianarisF(kvaziNiuntono), 13)}".ljust(20) +
          f" | {epsilonG}".ljust(10) +
          f" | {iterationKvaz}".ljust(15) + " |")
print('-' * 107)

# 1.4
# Gautas šaknų reikšmes patikrinkite naudodami išorinius išteklius
# ir pateikite patikrinimo rezultatus

print("\n1.4")
# daugianario saknis naudojant numpy.roots()
saknysNumpyF = numpy.roots(koeficientai)
print(f"Daugianario saknys su numpy.roots(): {saknysNumpyF}")

# funkcijos saknys naudojant scipy.optimize.fsolve()
saknysFsolveG = []
for saknuPora in funkcijosSaknys:
    saknisFsolve = fsolve(funkcijaG, (saknuPora[0] + saknuPora[1]) / 2)[0]
    saknysFsolveG.append(saknisFsolve)
print(f'Funkcijos saknys su scipy.optimize.fsolve: {saknysFsolveG}')

# 2.1
# Pateikite tarpinius grafikus, kai drauge su pateikta funkcija ℎ(𝑥 )
# nurodytame intervale atvaizduojama TE, kaijos narių skaičius lygus 3, 4 ir 5
