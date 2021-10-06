
import numpy as np


L = 20
dx = 0.2

def fD(x, x0, l):
    X = (x-x0)/l
    return np.cosh(X)**(-2)

def fB(x, x0, l):
    X = (x-x0)/l
    return 0.5*(1.+np.tanh(X))

def density(x):
    return 0.5+fD(x, 0.25*L, 1.0)+fD(x, 0.75*L, 1.0)

def by(x):
    return -1+2*(fB(x, 0.25*L, 1.0) - fB(x, 0.75*L, 1.0))

def b2(x):
    return by(x)**2

def T(x):
    K = 1
    return 1./density(x)*(K - 0.5*b2(x))

def vthx(x):
    return T(x)

