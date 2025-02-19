import numpy
import awkward

def deltaphi(one, another):
    a = one.phi
    b = another.phi
    return (a - b + numpy.pi) % (2 * numpy.pi) - numpy.pi

def delta_r(one, another):

    return numpy.hypot(one.eta - another.eta, deltaphi(one, another))



