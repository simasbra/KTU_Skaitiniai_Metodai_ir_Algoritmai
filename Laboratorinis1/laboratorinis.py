import sympy
import numpy
import math
import matplotlib.pyplot as pyplot
from scipy.optimize import fsolve
from sympy.plotting import plot


# -0.45x^4 + 1.04x^3 + 1.42x^2 - 2.67x - 0.97
def daugianaris_f(x):
    return -0.45 * x**4 + 1.04 * x**3 + 1.42 * x**2 - 2.67 * x - 0.97


# cos(2x) / (sin(x) + 1.5)) - (x / 5), kai -5 <= x <= 5
def funkcija_g(x):
    return (numpy.cos(2 * x) / (numpy.sin(x) + 1.5)) - (x / 5)


# -61cos(x) + cos(2x) + 12, kai 1 <= x <= 6
def funkcija_h(x):
    return -61 * numpy.cos(x) + numpy.cos(2 * x) + 12


# grubaus intervalo radimo funkcija
def grubus_intervalas(koeficientai):
    max_koef = max(koeficientai[:-1], key=abs)
    R = 1 + (max_koef / koeficientai[-1])
    return round(R, 2)


# tikslesnio intervalo radimo funkcija
def tikslesnis_intervalas(koeficientai):
    n = len(koeficientai)
    max_abs_index, B = min(enumerate(koeficientai[:-1]), key=lambda x: x[1])
    k = n - max_abs_index
    R = 1 + (abs(B) / koeficientai[-1])**(1 / k)
    return round(R, 2)


# skenavimo algoritmas
def skenavimas(funkcija, intervalas, zingsnis):
    # pradinė x reikšmė - intervalo pradžia
    x = intervalas[0]
    # čia bus laikomi šaknų atskyrimo intervalai
    saknys = []
    # nusprendžiama, kurią funkciją naudosime (pagal funkcijos parametrą)
    func = daugianaris_f if funkcija == 'f' else funkcija_g
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


# pusiaukirtos metodas
def pusiaukirtos(funkcija, saknys, tolerancija, iteracijosMax=100):
    patikslintosSaknys = []
    iterations = []
    for saknuPora in saknys:
        xFrom = saknuPora[0]
        xTo = saknuPora[1]
        xMid = (xFrom + xTo) / 2
        iteration = 0
        while (numpy.abs(funkcija(xMid)) > tolerancija):
            if (numpy.sign(funkcija(xMid)) == numpy.sign(funkcija(xFrom))):
                xFrom = xMid
            else:
                xTo = xMid
            xMid = (xFrom + xTo) / 2
            iteration += 1
        patikslintosSaknys.append(xMid)
        iterations.append(iteration)
        # print(f'Artinys ({xFrom}, {xTo})')
        # print(f'Patikslinta saknis: {xMid} Iteracijos: {format(iteration)}')
    return patikslintosSaknys, iterations


# kvazi niutono kirstiniu metodas
def kvazi_niutono(funkcija, saknys, tolerancija):
    patikslintosSaknys = []
    iterations = []
    for saknuPora in saknys:
        x0 = saknuPora[0]
        x1 = saknuPora[1]
        x2 = x1 - funkcija(x1) * (x1 - x0) / \
            float((funkcija(x1) - funkcija(x0)))
        iteration = 0
        while abs(funkcija(x2)) > tolerancija:
            funkcijax0 = funkcija(x0)
            funkcijax1 = funkcija(x1)
            if funkcijax1 == funkcijax0:  # Patikrina, ar yra dalyba iš nulio
                break
            x0, x1 = x1, x2
            iteration += 1
        patikslintosSaknys.append(x2)
        iterations.append(iteration)
    return patikslintosSaknys, iterations


# tikrina ar saknys skiriasi daugiau nei 1e-4
def ar_saknys_tinka(saknysH, saknysTE):
    for i in range(min(len(saknysH), len(saknysTE))):
        if abs(saknysH[i][0] - saknysTE[i][0]) > 1e-4:
            return False
    return True


def sudaryti_grafika_TE(f, fpirmas, x, i, saknys):
    numF = sympy.lambdify(x, f, "numpy")
    numFpirmas = sympy.lambdify(x, fpirmas, "numpy")
    valuesX = numpy.linspace(1, 6, 1000)
    valuesF = numF(valuesX)
    valuesFpirmas = numFpirmas(valuesX)

    pyplot.plot(valuesX, valuesF, color='r', linewidth=2)
    pyplot.plot(valuesX, valuesFpirmas, color='b', linewidth=0.5)
    pyplot.title(f"TE narių skaičius: {i}")
    pyplot.xlim([0.5, 6.5])
    pyplot.axhline(y=0, c='k', lw=1)

    for saknuPora in saknys:
        pyplot.plot(saknuPora[0], saknuPora[1], "og")
    pyplot.show()


# Funkcija skirtumams apskaiciuoti
def apskaiciuoti_skirtuma(saknysH, saknysTE):
    skirtumai = []
    for i in range(min(len(saknysH), len(saknysTE))):
        skirtumas = abs(saknysH[i] - saknysTE[i])
        skirtumai.append(skirtumas)
    return skirtumai


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
figureF = daugianaris_f(iverciaiF)
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
figureG = funkcija_g(intervalaiG)
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
epsilonF1 = 1e-8
pusiaukirtosF, iterationsPusF = pusiaukirtos(
    daugianaris_f, daugianarioSaknys, epsilonF1, 100)

epsilonG1 = 1e-8
pusiaukirtosG, iterationsPusG = pusiaukirtos(
    funkcija_g, funkcijosSaknys, epsilonG1, 100)

# kvazi-niutono (kirstiniu) metodas
epsilonF2 = 1e-13
kvaziNiutonoF, iterationsKvazF = kvazi_niutono(
    daugianaris_f, daugianarioSaknys, epsilonF2)

epsilonG2 = 1e-13
kvaziNiutonoG, iterationsKvazG = kvazi_niutono(
    funkcija_g, funkcijosSaknys, epsilonG2)

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
    kvaziNiutono = kvaziNiutonoF[i]
    iterationKvaz = iterationsKvazF[i]
    print("| Pusiaukirtos".ljust(15) +
          f" | [{round(saknuPora[0], 3)} - {round(saknuPora[1], 3)}]".ljust(25) +
          f" | {round(pusiaukirtos, 10)}".ljust(20) +
          f" | {round(daugianaris_f(pusiaukirtos), 5)}".ljust(20) +
          f" | {epsilonF1}".ljust(10) +
          f" | {iterationPus}".ljust(15) + " |")
    print("| Kvazi-Niutono".ljust(15) +
          f" | [{round(saknuPora[0], 3)} - {round(saknuPora[1], 3)}]".ljust(25) +
          f" | {round(kvaziNiutono, 10)}".ljust(20) +
          f" | {round(daugianaris_f(kvaziNiutono), 5)}".ljust(20) +
          f" | {epsilonF2}".ljust(10) +
          f" | {iterationKvaz}".ljust(15) + " |")
print('-' * 107)
print("| g(x) = cos(2x) / (sin(x) + 1.5)) - (x / 5), kai -5 <= x <= 5" +
      ' ' * 44 + '|')
print('-' * 107)
for i in range(len(funkcijosSaknys)):
    saknuPora = funkcijosSaknys[i]
    pusiaukirtos = pusiaukirtosG[i]
    iterationPus = iterationsPusG[i]
    kvaziNiutono = kvaziNiutonoG[i]
    iterationKvaz = iterationsKvazG[i]
    print("| Pusiaukirtos".ljust(15) +
          f" | [{round(saknuPora[0], 3)} - {round(saknuPora[1], 3)}]".ljust(25) +
          f" | {round(pusiaukirtos, 10)}".ljust(20) +
          f" | {round(funkcija_g(pusiaukirtos), 5)}".ljust(20) +
          f" | {epsilonG1}".ljust(10) +
          f" | {iterationPus}".ljust(15) + " |")
    print("| Kvazi-Niutono".ljust(15) +
          f" | [{round(saknuPora[0], 3)} - {round(saknuPora[1], 3)}]".ljust(25) +
          f" | {round(kvaziNiutono, 10)}".ljust(20) +
          f" | {round(funkcija_g(kvaziNiutono), 5)}".ljust(20) +
          f" | {epsilonG2}".ljust(10) +
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
    saknisFsolve = fsolve(funkcija_g, (saknuPora[0] + saknuPora[1]) / 2)[0]
    saknysFsolveG.append(saknisFsolve)
print(f'Funkcijos saknys su scipy.optimize.fsolve: {saknysFsolveG}')

# 2
# pateiktą funkciją ℎ(𝑥) išskleiskite Teiloro eilute (TE) nurodyto intervalo
# vidurio taško aplinkoje. Nustatykite TE narių skaičių, su kuriuo visos
# TE šaknys esančios nurodytame intervale, skiriasi nuo funkcijos ℎ(𝑥) šaknų
# ne daugiau negu |1e-4|. Tiek pateiktos funkcijos ℎ(𝑥) šaknis, tiek TE šaknist
# raskite antru iš pirmoje dalyje realizuotų skaitinių metodų
# (Niutono arba Kvazi-Niutono, priklausomai nuo varianto).aško aplinkoje.

# Nustatykite TE narių skaičių, su kuriuo visos TE šaknys esančios nurodytame intervale,
# nurodytame intervale atvaizduojama TE, kaijos narių skaičius lygus 3, 4 ir 5
# -61cos(x) + cos(2x) + 12, kai 1 <= x <= 6
print("2")
x, f, funk, fpirmas, df = sympy.symbols(('x', 'f', 'funk', 'fp', 'df'))
intervalasH = numpy.linspace(1, 6, 1000)
# funkcijos israika
f = -61 * sympy.cos(x) + sympy.cos(2 * x) + 12
funk = -61 * sympy.cos(x) + sympy.cos(2 * x) + 12
# pradinis artinys, TE pradinis taskas
x0 = (intervalasH[1] + intervalasH[-1]) / 2
# pirmas TE narys
fpirmas = f.subs(x, x0)

# funkcijos h grafikas
plotH1 = plot(f, (x, 1, 6), line_color='r', axis_center=(3, 4),
              show=False, title="TE narių skaičius: 0")
# pirmojo te nario grafikas
plotH2 = plot(fpirmas, (x, 1, 6), line_color='b',
              axis_center=(3, 4), show=False)
plotH1.append(plotH2[0])
plotH1.show()

# h(x) funkcijos saknu atskyrimo intervalai
saknysH, iteracijosH = skenavimas(funkcija_h, (1, 6), 0.01)
# h(x) funkcijos saknys kvazi niutono kirstiniu metodu
kvaziNiutonoH = kvazi_niutono(funkcija_h, saknysH, 1e-10)
saknuSkaicius = [(0, 0)]
skirtumaiPerIteracija = [[] for _ in range(len(kvaziNiutonoH[0]))]
# iki kokios eiles naudoti TE narius
N = 100
for i in range(1, N + 1):
    # analizinis diferencijavimas
    funk = funk.diff(x)
    # TE suma
    fpirmas = fpirmas + funk.subs(x, x0) / math.factorial(i) * (x - x0)**i
    saknysTE, iteracijosTE = skenavimas(
        sympy.lambdify(x, fpirmas, "numpy"), (1, 6), 0.01)
    kvaziNiutonoTE = kvazi_niutono(
        sympy.lambdify(x, fpirmas, "numpy"), saknysTE, 1e-10)
    saknuSkaicius.append((i, len(kvaziNiutonoH)))
    skirtumai = apskaiciuoti_skirtuma(kvaziNiutonoH[0], kvaziNiutonoTE[0])
    for j in range(len(skirtumai)):
        skirtumaiPerIteracija[j].append(skirtumai[j])
    if i >= 3 and i <= 5:
        print(f"Te nariu: {i}, funkcija {fpirmas}")
        sudaryti_grafika_TE(f, fpirmas, x, i, kvaziNiutonoH)
    if ar_saknys_tinka(kvaziNiutonoH, kvaziNiutonoTE):
        print(f"\nvisos saknu reiksmes intervale skiriasi maziau nei {
              1e-4}. Is viso yra {i} nariu")
        sudaryti_grafika_TE(f, fpirmas, x, i, kvaziNiutonoH)
        break

# daugianaris simboliais
a = sympy.Poly(fpirmas, x)
print(f"\ndaugianaris simboliais {a}")
# visi koeficientai nuo vyriausio
kf = numpy.array(a.all_coeffs())
saknys = numpy.roots(kf)
print(f"\nvisi koeficientai nuo vyriausio {saknys}")

# Grafikas, rodantis saknu skaiciu priklausomai nuo TE nariu skaiciaus
pyplot.figure(figsize=(10, 6))
pyplot.plot([x[0] for x in saknuSkaicius], [x[1] for x in saknuSkaicius])
pyplot.title('Šaknų skaičius priklausomai nuo TE narių skaičiaus')
pyplot.xlabel('TE narių skaičius')
pyplot.ylabel('Šaknų skaičius')
pyplot.show()

# atskiri skirtumo grafikai kiekvienai sakniai
for idx, skirtumai in enumerate(skirtumaiPerIteracija):
    pyplot.figure(figsize=(10, 6))
    pyplot.plot(range(1, len(skirtumai) + 1), skirtumai, marker='o')
    pyplot.title(f'Šaknies {idx + 1} tikslumo skirtumai tarp h(x) ir TE')
    pyplot.xlabel('TE narių skaičius')
    pyplot.ylabel('Skirtumas tarp h(x) ir TE šaknies')
    pyplot.grid(True)
    pyplot.show()
