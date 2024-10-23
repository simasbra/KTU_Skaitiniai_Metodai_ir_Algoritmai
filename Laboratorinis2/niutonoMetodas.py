import numpy


def LF(x):
    s = numpy.array([x[0]**2 + (x[1] + numpy.cos(x[0]))**2 - 40,
                     (x[0]/2)**3 + 25*x[1]**2 - 50])
    return s


def nuliai(x1, x2):
    result = numpy.array([[x1**2 + (x2 + numpy.cos(x1))**2 - 40],
                          [(x1 / 2)**3 + 25*x2**2 - 50]])
    return result


def dZ(x1, x2):
    # dZ nuo x1, dZ nuo x2
    result = numpy.array([
        [2 * x1 - 2 * (numpy.cos(x1) + x2) * numpy.sin(x1), 2 * (x2 + numpy.cos(x1))],
        [3 * x1**2 / 8, 50 * x2]])
    return result


def niutono_metodas(x1, x2, iteracijos):
    x = -1 * numpy.ones((2, 1))
    x[0] = x1
    x[1] = x2

    for i in range(iteracijos):
        try:
            x = x - numpy.matmul(numpy.linalg.inv(
                dZ(x[0, 0], x[1, 0])), nuliai(x[0, 0], x[1, 0]))
            Zx = numpy.linalg.norm(nuliai(x[0, 0], x[1, 0]))
        except Exception:
            # print(f"Matrica singuliari kai x1 = {x[0]} x2 = {x[1]}")
            return (x1, x2, None, None)
        # print(f"\nIteracija {i}: f(x) = {Zx}")

        if abs(Zx) < 1e-10:
            break
    # print(f"RASTA GALUTINE REIKSME: {Zx}\n x1 = {x[0]}, x2 = {x[1]}")
    # print(f"Iteraciju skaicius: {i}")
    return (x1, x2, x[0][0], x[1][0])


def niutono_sprendiniai(grid_x, grid_y):
    sprendiniai = []
    maxIteracijos = 100
    for i in range(len(grid_x)):
        for j in range(len(grid_y)):
            spr = niutono_metodas(grid_x[i, j], grid_y[i, j], maxIteracijos)
            sprendiniai.append(spr)
    return sprendiniai


# Priskirti spalvas kiekvienam sprendiniui pagal jo artuma
# filtruotam sprendiniui
def niutono_spalvos(sprendiniai, filtSprendiniai, spalvos, eps, defaultSpalva):
    priskirtosSpalvos = []
    for spr in sprendiniai:
        if spr[2] is None:
            priskirtosSpalvos.append(defaultSpalva)
            continue
        # Spalvos pagal artuma prie filtruotu sprendiniu
        for idx, egzSprendiniai in enumerate(filtSprendiniai):
            if abs(egzSprendiniai[2] - spr[2]) < eps and abs(
                    egzSprendiniai[3] - spr[3]) < eps:
                # Iteruoti per spalvas
                priskirtosSpalvos.append(spalvos[idx % len(spalvos)])
                break
    return priskirtosSpalvos


# Filtruoti sprendinius pagal jo artuma su paklaida
def filtruoti_sprendinius(sprendiniai, eps):
    filtSprendiniai = []
    for sprendinys in sprendiniai:
        if sprendinys[2] is None:
            continue
        if not any(abs(existing_spr[2] - sprendinys[2]) < eps and abs(
            existing_spr[3] - sprendinys[3]) < eps
                for existing_spr in filtSprendiniai):
            filtSprendiniai.append(sprendinys)
    return filtSprendiniai
