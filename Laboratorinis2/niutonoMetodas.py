import numpy


def funkcija_Z1(x):
    return x[0]**2 + (x[1] + numpy.cos(x[0]))**2 - 40


def funkcija_Z2(x):
    return (x[0]/2)**3 + 25*x[1]**2 - 50


def netiesiniu_lygciu_sistema(x):
    return numpy.array([funkcija_Z1(x), funkcija_Z2(x)])
