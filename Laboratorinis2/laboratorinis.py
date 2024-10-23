import numpy


def atspindys_su_qr(A, b):
    # lygciu skaicius nustatomas pagal ivesta matrica A
    n = (numpy.shape(A))[0]
    # laisvuju nariu vektoriu skaicius nustatomas pagal ivesta matrica b
    nb = (numpy.shape(b))[1]
    print("fu")


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

atspindys_su_qr(A7, b7)

