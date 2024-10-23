import numpy


def atspindys_su_qr(A, b):
    # lygciu skaicius nustatomas pagal ivesta matrica A
    n = (numpy.shape(A))[0]
    # laisvuju nariu vektoriu skaicius nustatomas pagal ivesta matrica b
    nb = (numpy.shape(b))[1]
    # isplestoji matrica
    A1 = numpy.hstack((A, b))

    det = numpy.linalg.det(A)
    if det == 0:
        print("Matrica yra singuliari")

    # tiesioginis etapas(atspindziai):
    for i in range(0, n - 1):
        z = A1[i:n, i]
        zp = numpy.zeros(numpy.shape(z))
        zp[0] = numpy.linalg.norm(z)
        omega = z - zp
        omega = omega/numpy.linalg.norm(omega)
        Q = numpy.identity(n - i) - 2 * omega * omega.transpose()
        A1[i:n, :] = Q.dot(A1[i:n, :])
    print(f"Gauta matrica: \n {A1.round(2)}")

    # atgalinis etapas:
    x = numpy.zeros(shape=(n, nb))
    eps = 1e-15

    if abs(A1[n - 1, n]) < eps and abs(A1[n - 1, n - 1]) < eps:
        x[n - 1] = 1
        print(f"Kintamasis x{n} gali buti bet koks skaicius")
        print(f"Tarkime, x{n} = 1")
    elif abs(A1[n - 1, n]) > eps and abs(A1[n - 1, n - 1]) < eps:
        print("Sprendiniu nera")
        return
    else:
        x[n - 1] = A1[n - 1, n] / A1[n - 1, n - 1]

    # range pradeda n-1 ir baigia 0 (trecias parametras yra zingsnis)
    for i in range(n - 2, -1, -1):
        r = A1[i, n] - A1[i, i+1:n] * x[i+1:n]
        if A1[i, i] == 0 and abs(r) < eps:
            x[i] = 1
            print(f"r = {r}\nKintamasis x{i+1} gali buti bet koks skaicius")
            print(f"Tarkime, x[{i+1}] = 1:")
        elif A1[i, i] == 0 and abs(r) > eps:
            print(f"r = {r}\nKintamasis x{i+1}")
            print("Sprendiniu nera")
            return
        else:
            x[i] = r / A1[i, i]
            print(f"x{i+1} = {x[i][0]}")

    print("Gautas rezultatas:")
    print(f"x =\n{x.round(2)}")
    print(f"A =\n{A.round(2)}")
    print(f"Sprendinio patikrinimas:\n{(A.dot(x) - b).round(2)}")

    print("Sprendinys pagal python biblioteka:")
    try:
        print(numpy.linalg.solve(A, b))
    except Exception:
        print("Pagal python biblioteka nustatyta, kad matrica yra singuliari")


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
