# Importuojame reikalingas bibliotekas
import sympy as sp  # simbolinės matematikos biblioteka
import numpy as np   # skaičiavimams ir duomenų manipuliacijai
import math          # matematinė funkcijos
import matplotlib.pyplot as plt  # grafikų braižymui

# Apibrėžiame funkciją, kuri grąžina sinusą


def funk(x): return np.sin(x)


# Nustatome Pi vertę
Pi = math.pi

# Sukuriame x vertes intervale nuo -2π iki 2π, su 50 taškų
xxx = np.linspace(-2*Pi, 2*Pi, 50)
# Apskaičiuojame atitinkamas y vertes naudodami funk(x)
fff = funk(xxx)
# Atspausdiname y vertes
print(fff)

# Braižome funkcijos grafiką
plt.plot(xxx, fff, 'b-')  # 'b-' - mėlyna linija
plt.grid()  # Pridėti tinklelio linijas
plt.xlabel('x')  # X ašies pavadinimas
plt.ylabel('f')  # Y ašies pavadinimas

# Sukuriame naują figūrą ir dvi subfigūras
fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)  # Pirma subfigūra
ax2 = fig.add_subplot(1, 2, 2)  # Antra subfigūra
f1, = ax1.plot(xxx, fff, 'b-.')  # Mėlyna linija su punktyru
ax1.grid()  # Tinklelis pirmoje subfigūroje
ax1.set_xlabel('x')  # X ašies pavadinimas
ax1.set_ylabel('f')  # Y ašies pavadinimas

# Antroje subfigūroje braižome f(x) - x
f2, = ax2.plot(xxx, fff - xxx, 'r-*')  # Raudoni žvaigždutės taškai
ax2.grid()  # Tinklelis antroje subfigūroje
ax2.set_xlabel('x')  # X ašies pavadinimas
ax2.set_ylabel('f')  # Y ašies pavadinimas

# Simbolinė išraiška duotos funkcijos grafikui
x, f = sp.symbols(('x', 'f'))  # Nustatome simbolius x ir f
print(x, f)  # Atspausdiname simbolius
f = sp.sin((x**2))  # Apibrėžiame funkciją f kaip sin(x^2)
print(f)  # Atspausdiname funkciją

# Sukuriame x vertes intervale nuo -2π iki 2π su 500 taškų
xxx = np.linspace(-2*Pi, 2*Pi, 500)
# Inicializuojame tuščią masyvą y vertėms
fff = np.array(np.zeros(np.size(xxx)))
# Apskaičiuojame y vertes naudojant simbolinę funkciją
for i in range(len(xxx)):
    fff[i] = f.subs(x, xxx[i]).evalf()  # Apskaičiuojame f(x)

# Braižome simbolinės funkcijos grafiką
plt.plot(xxx, fff, 'b-')  # Mėlyna linija
plt.grid()  # Tinklelis
plt.xlabel('x')  # X ašies pavadinimas
plt.ylabel('f')  # Y ašies pavadinimas
plt.show()  # Rodome grafiką

# Apskaičiuojame funkcijos išvestinę
df = sp.diff(f, x)  # Apskaičiuojame išvestinę
# Inicializuojame tuščią masyvą išvestinėms
dfff = np.array(np.zeros(np.size(xxx)))
# Apskaičiuojame išvestines vertes
for i in range(len(xxx)):
    dfff[i] = df.subs(x, xxx[i]).evalf()  # Apskaičiuojame df/dx

# Braižome išvestinės grafiką
plt.plot(xxx, dfff, 'r-')  # Raudona linija
plt.grid()  # Tinklelis
plt.xlabel('x')  # X ašies pavadinimas
plt.ylabel('df')  # Y ašies pavadinimas

# Šaknų atskyrimo metodas: pusiaukirta

# Apibrėžiame funkciją, kuri grąžina sinusą (dar kartą)


def funk(x): return np.sin(x)


# Nustatome pradinius intervalus
Pi = math.pi
xn = -6  # Kairysis intervalas
xn1 = 2  # Dešinysis intervalas

# Sukuriame x vertes intervale nuo -2π iki 2π su 50 taškų
xxx = np.linspace(-2*Pi, 2*Pi, 50)
# Apskaičiuojame y vertes
fff = funk(xxx)

# Nustatome pradinį tikslumą
tikslumas = 9999
it = 0  # Iteracijų skaičius
# Pusiaukirtos metodas
while tikslumas > 1e-6:  # Kol tikslumas didesnis nei 1e-6
    it += 1  # Padidiname iteracijų skaičių
    fn = funk(xn)  # Apskaičiuojame funkcijos vertę kairiajame intervale
    fn1 = funk(xn1)  # Apskaičiuojame funkcijos vertę dešiniajame intervale
    xmid = (xn + xn1) / 2  # Apskaičiuojame vidurį
    fmid = funk(xmid)  # Apskaičiuojame funkcijos vertę viduryje

    # Braižome grafiką su šaknimis
    plt.plot(xxx, fff, 'b-.')  # Mėlyna linija
    plt.grid()  # Tinklelis
    plt.xlabel('x')  # X ašies pavadinimas
    plt.ylabel('f')  # Y ašies pavadinimas
    plt.plot(xn, 0, 'k*')  # Juoda žvaigždutė kairiajame intervale
    plt.plot(xn1, 0, 'b*')  # Mėlyna žvaigždutė dešiniajame intervale
    plt.plot(xmid, 0, 'rp')  # Raudona punktyrinė žvaigždutė viduryje
    # Juoda brūkšninė linija nuo kairiojo intervalo
    plt.plot([xn, xn], [0, fn], 'k--')
    # Juoda brūkšninė linija nuo dešiniojo intervalo
    plt.plot([xn1, xn1], [0, fn1], 'k--')

    # Tikriname, kurioje pusėje yra šaknis
    if np.sign(fn) == np.sign(fmid):
        xn = xmid  # Jei kairėje pusėje, atnaujiname kairįjį intervalą
    else:
        xn1 = xmid  # Jei dešinėje pusėje, atnaujiname dešinįjį intervalą
    tikslumas = np.abs(xn - xn1)  # Apskaičiuojame tikslumą
    # Atspausdiname informaciją apie iteracijas
    print('iteracija ', it, 'saknis= ', xmid, ' tikslumas= ', tikslumas)
    plt.show()  # Rodome grafiką
