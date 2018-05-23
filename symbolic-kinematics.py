from sympy import *
import numpy

class c(Function):

    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x is S.Zero:
                return cos(x)

class s(Function):

    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x is S.Zero:
                return sin(x)

class b(Function): #alfa
    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x is 2: return  sympy.pi/2
            elif x is 4: return -sympy.pi/2
            elif x is 5: return sympy.pi/2
            elif x is 6: return -sympy.pi/2
            else: return 0

class a(Function):
    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x is 1: return 'a1'
            if x is 2: return 'a2'
            if x is 3: return 'a3'
            else: return 0

class d(Function):
    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x is 4: return 'd4'
            if x is 7: return 'd7'
            else: return 0

def getmatrix(i):
    mat = Matrix([
        [c(i), -s(i), 0, a(i-1)],
        [s(i)*c(b(i-1)), c(i)*c(b(i-1)), -s(b(i-1)), -d(i)*s(b(i-1))],
        [s(i)*s(b(i-1)), c(i)*s(b(i-1)), c(b(i-1)), d(i)*c(b(i-1))],
        [0,0,0,1]
    ])
    return mat;
    


init_printing(use_unicode=True)
T1 = getmatrix(1)
T2 = getmatrix(2)
T3 = getmatrix(3)
T4 = getmatrix(4)
T5 = getmatrix(5)
T6 = getmatrix(6)
pprint(simplify(T1*T2))
