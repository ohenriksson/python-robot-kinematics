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
        else:
            return cos(x)

class s(Function):
    @classmethod
    def eval(cls,x):
        if x == pi/2: return 1
        if x == -pi/2: return -1
        if x == 0:    return 0
        else:
            return sin(x)

class b(Function): #alfa
    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x == 1: return  -pi/2
            elif x == 3: return -pi/2
            elif x == 4: return pi/2
            elif x == 5: return -pi/2
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
    return mat
    
def MyInv(matrix):
    myinv = matrix
    A = matrix[0:3,0:3]
    R = matrix[0:3,3]
    Rnew = -(A.T)*R
    myinv[0:3,0:3] = simplify(A.T)
    myinv[0:3,3] = Rnew
    return myinv 



def printall():
    for x in range(1,7):
        str1 = '\\begin{equation} \prescript{' + str(x-1) +'}{' +str(x) +'}{DH} = '
        str2 = '\end{equation}'
        # print(str1 + latex(getmatrix(x)) +str2)
        pprint(getmatrix(x))

def createPosMatrix():
    x, y, z = symbols('x y z')
    Pmat = eye(4)
    Pmat[0:3,3] = Matrix([[x,y,z]]).T
    return Pmat

def solveTheta1(Pmat2, T16):
    print('-----theta1-------')
    py = T16.row(1).col(3)[0]
    temp1 = Pmat2.row(1).col(3)[0] - py
    theta1 = solve(temp1,t(1))
    pprint(theta1)

def solveTheta3(Pmat2, T16):
    print('-----theta3-------')
    px = T16.row(0).col(3)[0]
    py = T16.row(1).col(3)[0]
    pz = T16.row(2).col(3)[0]
    eq1 =  - pow(px,2) + pow(Pmat2.row(0).col(3)[0],2)
    eq2 =  - pow(py,2) + pow(Pmat2.row(1).col(3)[0],2)
    eq3 =  - pow(pz,2) + pow(Pmat2.row(2).col(3)[0],2)
    pprint(simplify(eq1))
    pprint(simplify(eq2))
    pprint(simplify(eq3))
    print('-----')
    slask = simplify(-eq1-eq2-eq3)
    pprint(slask)
    theta3 = solve(slask,t(3))
    pprint(theta3)

    # pprint(simplify(eq1 - eq2 -eq3 ))

def solveTheta2(T03, T06, T36):
    leftM = MyInv(T03)*createPosMatrix()
    l14 = leftM.row(0).col(3)[0]
    r14 = T36.row(0).col(3)[0]
    l24 = leftM.row(1).col(3)[0]
    r24 = T36.row(1).col(3)[0]
    eq1 = l14-r14
    eq2 = l24-r24
    sin23 = solve(eq1,sin(t(2)+t(3)))
    pprint(sin23)
    cos23 = solve(eq2,cos(t(2)+t(3)))
    pprint(cos23)





#matrices
T1 = getmatrix(1)
T2 = getmatrix(2)
T3 = getmatrix(3)
T4 = getmatrix(4)
T5 = getmatrix(5)
T6 = getmatrix(6)

init_printing(use_unicode=True)

#the T16 matrix:

# #the real T16 matrix
# T16 = simplify(T2*T3*T4*T5*T6)

    # Pmat2 = MyInv(T1)*Pmat
# solveTheta1(Pmat2, T16)
# solveTheta3(Pmat2, T16)
T03 = T1*T2*T3
T06 = T1*T2*T3*T4*T5*T6
T36 = T4*T5*T6
solveTheta2(T03,T06,T36)

# T16a = T1.inv()*(T1*T2*T3*T4*T5*T6)
# T16b = T2*T3*T4*T5*T6
# cella =(T16a.col(3).row(1))
# cellb =(T16b.col(3).row(1))


# pprint(simplify(cellb))
