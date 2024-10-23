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
