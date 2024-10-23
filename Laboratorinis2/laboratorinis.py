import numpy
from atspindysSuQR import atspindys_su_qr

# 1 Tiesinių lygčių sistemų sprendimas

# Lentelėje 1 duotos tiesinės lygčių sistemos, 2 lentelėje nurodytas metodas ir
# lygčių sistemų numeriai (iš 1 lentelės). Reikia suprogramuoti nurodytą metodą
# ir išspręsti pateiktas lygčių sistemas. Programoje sprendžiant lygčių
# sistemas turi būti įvertinti atvejai:
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
# c. Nagrinėjamoje srityje sudarykite stačiakampį tinklelį (𝑥1, 𝑥2 poras).
# Naudodami užduotyje nurodytą metodą apskaičiuokite netiesinių lygčių sistemos
# sprendinius, kai pradinis artinys įgyja tinklelio koordinačių reikšmes.
# Tinklelyje vienodai pažymėkite taškus, kuriuos naudojant kaip pradinius
# artinius gaunamas tas pats sprendinys. Lentelėje pateikite apskaičiuotus skirtingus
# sistemos prendinius ir bent po vieną jam atitinkantį pradinį artinį.
# d. Gautus sprendinius patikrinkite naudodami išorinius išteklius

