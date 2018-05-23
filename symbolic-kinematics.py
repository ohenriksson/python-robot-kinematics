from sympy import *
from sympy.abc import *
import mpmath
import numpy

mpmath.mp.pretty = True

class c(Function):
    @classmethod
    def eval(cls,x):
        if x == pi/2: return 0
        if x == -pi/2: return 0
        if x == 0:    return 1
        if x.is_Number:
            return cos(x)

class s(Function):
    @classmethod
    def eval(cls,x):
        if x == pi/2: return 1
        if x == -pi/2: return -1
        if x == 0:    return 0
        if x.is_Number:
            return sin(x)

class b(Function): #alfa
    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x == 2: return  -pi/2
            elif x == 4: return -pi/2
            elif x == 5: return pi/2
            elif x == 6: return -pi/2
            else: return 0

class a(Function):
    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x == 1:   return Symbol('a1')
            elif x == 2: return Symbol('a2')
            elif x == 3: return Symbol('a3')
            else: return 0

class d(Function):
    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x == 4: return Symbol('d4')
            elif x == 7: return Symbol('d7')
            else: return 0

class t(Function):
    @classmethod
    def eval(cls,x):
        if x.is_Number:
            return Symbol(str(theta) +str(x))

def getmatrix(i):
    mat = Matrix([
        [c(t(i)), -s(t(i)), 0, a(i-1)],
        [s(t(i))*c(b(i-1)), c(t(i))*c(b(i-1)), -s(b(i-1)), -d(i)*s(b(i-1))],
        [s(t(i))*s(b(i-1)), c(t(i))*s(b(i-1)), c(b(i-1)), d(i)*c(b(i-1))],
        [0,0,0,1]
    ])
    return mat;
    

def printall():
    for x in range(1,7):
        str1 = '\\begin{equation} \prescript{' + str(x-1) +'}{' +str(x) +'}{DH} = '
        str2 = '\end{equation}'
        print(str1 + latex(getmatrix(x)) +str2)


init_printing(use_unicode=True)
T1 = getmatrix(1)
T2 = getmatrix(2)
T3 = getmatrix(3)
T4 = getmatrix(4)
T5 = getmatrix(5)
T6 = getmatrix(6)

T16a = T1.inv()*(T1*T2*T3*T4*T5*T6)
T16b = T2*T3*T4*T5*T6
cella =(T16a.col(3).row(1))
cellb =(T16b.col(3).row(1))


pprint(simplify(cellb))
