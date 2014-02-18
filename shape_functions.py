import numpy as np
from numpy.linalg import *


def PHI0(xi,eta,mu):
    return 0.125*(1-xi)*(1-eta)*(1-mu) 

def PHI1(xi,eta,mu):
    return 0.125*(1+xi)*(1-eta)*(1-mu)

def PHI2(xi,eta,mu):
    return 0.125*(1+xi)*(1+eta)*(1-mu) 

def PHI3(xi,eta,mu):
    return 0.125*(1-xi)*(1+eta)*(1-mu) 

def PHI4(xi,eta,mu):
    return 0.125*(1-xi)*(1-eta)*(1+mu)

def PHI5(xi,eta,mu):
    return 0.125*(1+xi)*(1-eta)*(1+mu) 

def PHI6(xi,eta,mu):
    return 0.125*(1+xi)*(1+eta)*(1+mu) 

def PHI7(xi,eta,mu):
    return 0.125*(1-xi)*(1+eta)*(1+mu) 



def dPHI0(xi,eta,mu):
    return [ 0.125*(-1)*(1-eta)*(1-mu),\
    0.125*(1-xi)*(-1)*(1-mu),\
    0.125*(1-xi)*(1-eta)*(-1)]

def dPHI1(xi,eta,mu):
    return  [0.125*(1)*(1-eta)*(1-mu), \
    0.125*(1+xi)*(-1)*(1-mu),\
    0.125*(1+xi)*(1-eta)*(-1)]

def dPHI2(xi,eta,mu):
    return [ 0.125*(1)*(1+eta)*(1-mu),\
    0.125*(1+xi)*(1)*(1-mu),\
    0.125*(1+xi)*(1+eta)*(-1) ]

def dPHI3(xi,eta,mu):
    return [ 0.125*(-1)*(1+eta)*(1-mu),\
    0.125*(1-xi)*(1)*(1-mu),\
    0.125*(1-xi)*(1+eta)*(-1)] 

def dPHI4(xi,eta,mu):
    return [0.125*(-1)*(1-eta)*(1+mu),\
    0.125*(1-xi)*(-1)*(1+mu),\
    0.125*(1-xi)*(1-eta)*(1)]

def dPHI5(xi,eta,mu):
    return [ 0.125*(1)*(1-eta)*(1+mu),\
    0.125*(1+xi)*(-1)*(1+mu),\
    0.125*(1+xi)*(1-eta)*(1)]

def dPHI6(xi,eta,mu):
    return [ 0.125*(1)*(1+eta)*(1+mu),\
    0.125*(1+xi)*(1)*(1+mu),\
    0.125*(1+xi)*(1+eta)*(1)]

def dPHI7(xi,eta,mu):
    return [ 0.125*(-1)*(1+eta)*(1+mu),\
    0.125*(1-xi)*(1)*(1+mu),\
    0.125*(1-xi)*(1+eta)*(1)]


