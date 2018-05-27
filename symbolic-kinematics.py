from sympy import *
from sympy.abc import *
import mpmath
import numpy

mpmath.mp.pretty = True
init_printing(use_unicode=True)

class c(Function):
    @classmethod
    def eval(cls,x):
        # return cos(x)
        if x == pi/2: return 0
        if x == -pi/2: return 0
        if x == 0:    return 1
        else:
            return cos(x)

class s(Function):
    @classmethod
    def eval(cls,x):
        # return sin(x)
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
        if x.is_String:
            return Symbol('a'+x)

class d(Function):
    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x == 4: return Symbol('d4')
            elif x == 7: return Symbol('d7')
            else: return 0
        else: return Symbol('d_'+str(x))

class t(Function):
    @classmethod
    def eval(cls,x):
        if x.is_Number:
            if x == 0: return 0
            else: return Symbol(str(theta) +str(x))
        else: return 0

def getmatrix(i, dval=0):
    if dval == 0: dval = d(i)
    else: dval = d(dval)
    mat = Matrix([
        [c(t(i)), -s(t(i)), 0, a(i-1)],
        [s(t(i))*c(b(i-1)), c(t(i))*c(b(i-1)), -s(b(i-1)), -dval*s(b(i-1))],
        [s(t(i))*s(b(i-1)), c(t(i))*s(b(i-1)), c(b(i-1)), dval*c(b(i-1))],
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



def printall(list1,indices):
    for x,i in zip(list1,indices):
        str1 = '\\begin{equation} \prescript{' + i[0] +'}{' +i[1] +'}{DH} = '
        str2 = '\end{equation}'
        print(str1 + latex(x) +str2)
        pprint(x)

def createPosMatrix():
    x, y, z = symbols('p_x p_y p_z')
    Pmat = eye(4)
    Pmat[0:3,3] = Matrix([[x,y,z]]).T
    return Pmat

def solveTheta1(T01, T16):
    print('-----theta1-------')
    leftM = MyInv(T01)*createPosMatrix()
    py = T16.row(1).col(3)[0]
    temp1 = Eq(leftM.row(1).col(3)[0],py)
    theta1 = solve(temp1,t(1))
    print(latex(theta1[0]))
    print(latex(theta1[1]))

def solveTheta3(T01, T16):
    p_x, p_y, p_z = symbols('p_x p_y p_z')
    print('-----theta3-------')
    leftM = MyInv(T01)*createPosMatrix()
    px = T16.row(0).col(3)[0]
    py = T16.row(1).col(3)[0]
    pz = T16.row(2).col(3)[0]
    l14 = leftM.row(0).col(3)[0]
    l24 = leftM.row(1).col(3)[0]
    l34 = leftM.row(2).col(3)[0]

    #the solver could not solve it with a on the wrong side of the equal sign(probably because of the squaring later)
    l14 = l14 -a(1)
    px = px - a(1)

    eqx = simplify(Eq(l14,px)**2)
    eqy = simplify(Eq(l24,py)**2)
    eqz = simplify(Eq(l34,pz)**2)

    #sadly, squaring the whole Eq did not follow through, workaround:
    eqx = pow(l14,2)-pow(px,2)
    eqy = pow(l24,2)-pow(py,2)
    eqz = pow(l34,2)-pow(pz,2)

    eq4 = Eq(simplify(eqx+eqy+eqz),0)
    pprint(eq4)
    lhs = 2*a(2) * (a(3)*c(t(3)) - d(4)*s(t(3))) 
    rhs = a(1)**2 - a(3)**2 -d(4)**2 + p_x**2 + p_y**2 -a(2)**2 + p_z**2 - 2*a(1)*p_x*c(t(1)) - 2*a(1) *p_y *s(t(1)) 
    lhs = solve(eq4,rhs)[0]

    lhs = lhs/(2*a(2))
    K = rhs/(2*a(2))
    # slask = solve(eq4,     slask = solve(eq4, (d(4)*s(t(3))) )
    print(latex(lhs))
    print(latex(K))
    pprint(lhs)
    pprint(K)
    K = symbols('K')

    theta3a = atan2(a(3),d(4)) - atan2(K,+sqrt(a(3)**2+d(3)**2-K**2))
    theta3b = atan2(a(3),d(4)) - atan2(K,-sqrt(a(3)**2+d(3)**2-K**2))
    print(latex(theta3a))
    return theta3a


def solveTheta2(T03, T36):
    p_x, p_y, p_z = symbols('p_x p_y p_z')
    print('-----theta2-------')
    leftM = MyInv(T03)*createPosMatrix()
    l14 = leftM.row(0).col(3)[0]
    r14 = T36.row(0).col(3)[0]
    l24 = leftM.row(1).col(3)[0]
    r24 = T36.row(1).col(3)[0]
    eq1 = l14-r14
    eq2 = l24-r24

    denumA = c(t(1))*p_x + s(t(1))*p_y -a(1)
    denum = p_z**2 + denumA**2

    nom1 = p_z*(-a(3)-a(2)*c(t(3))) +  denumA*(a(2)* s(t(3))-d(4))
    eq1 = Eq(s(t(2)+t(3)),nom1/denum)
    nom2 = p_z*(a(2)*s(t(3))-d(4)) -  denumA*(a(3) + a(2)*c(t(3)))
    eq2 = Eq(c(t(2)+t(3)),nom2/denum)

    eq1 = (simplify(eq1))
    eq2 = (simplify(eq2))
    print(latex(eq1))
    print(latex(eq2))
    print('---')

    EQ = Eq(t(2)+t(3),atan2(nom1,nom2))
    print(latex(EQ))
    pprint(EQ)

    return
    # c1px = solve(eq1,c(t(1)*p_x)))[0]
    # c23 = solve(eq2,sin(t(2)+t(3)))[0]
    # pprint(simplify(s23))
    # pprint(simplify(c23))
    # print('---')

    # eq2b = eq2.subs(sin(t(2)+t(3)),s23)
    # pprint(simplify(eq2b))
    # return

    # A, B = fraction(solve(eq1,c(t(2)+t(3)))[0])
    # C, B = fraction(solve(eq2,s(t(2)+t(3)))[0])
    # #assume z>0
    # theta23 = atan2(A,C)
    # theta2 = theta23 - t(3)
    # pprint(theta23)
    # print('----')
    # pprint(latex(Eq(t(2)+t(3),theta23)))
    # print(latex(Eq(t(2)+t(3),theta23)))
    # return theta2


def solveTheta4(T03, T36):
    print('-----theta4-------')
    leftM = MyInv(T03)*createPosMatrix()
    c = 2
    l13 = leftM.row(0).col(c)[0]
    r13 = T36.row(0).col(c)[0]
    l33 = leftM.row(2).col(c)[0]
    r33 = T36.row(2).col(c)[0]
    eq1 = Eq(l13,r13)
    eq2 = Eq(l33,r33)

    print(latex(eq1))
    print(latex(eq2))

    c4 = solve(eq1, cos(t(4))*sin(t(5)))[0]
    s4 = solve(eq2, sin(t(4))*sin(t(5)))[0]

    #assume s5!=0 #singularity
    theta4 = Eq(t(4),atan2(s4,c4))
    pprint(theta4)
    print(latex(theta4))
    return theta4


def solveTheta5(T04, T46):
    print('-----theta5-------')
    leftM = MyInv(T04)*createPosMatrix()
    c = 2
    l13 = leftM.row(0).col(c)[0]
    r13 = T46.row(0).col(c)[0]
    l33 = leftM.row(2).col(c)[0]
    r33 = T46.row(2).col(c)[0]
    s5 = solve(l13-r13,sin(t(5)))[0]
    c5 = solve(l33-r33,cos(t(5)))[0]
    theta5 = Eq(t(5),atan2(s5,c5))

    pprint(theta5)
    print(latex(theta5))

    return theta5


def solveTheta6(T05, T56):
    print('-----theta6-------')
    leftM = MyInv(T05)*createPosMatrix()
    c = 0
    l11 = leftM.row(0).col(c)[0]
    r11 = T56.row(0).col(c)[0]
    l31 = leftM.row(2).col(c)[0]
    r31 = T56.row(2).col(c)[0]
    c6 = solve(l11-r11,cos(t(6)))[0]
    s6 = solve(l31-r31,sin(t(6)))[0]
    theta6 = Eq(t(6),atan2(s6,c6))

    pprint(theta6)
    print(latex(theta6))
    return theta6



#---static matrices
T0 = getmatrix(0,'b')
T7 = getmatrix(7)
# printall([T0, T7],[['b','0'],['6','7']])
# print('----')

#---dynamic matrices
T1 = getmatrix(1)
T2 = getmatrix(2)
T3 = getmatrix(3)
T4 = getmatrix(4)
T5 = getmatrix(5)
T6 = getmatrix(6)
# dhnumbers = map(lambda n: str(n), [[0,1],[1,2],[2,3],[3,4],[4,5],[5,6]])
# printall([T1, T2, T3, T4, T5, T6],dhnumbers)
# print('----')

#---inverse solving...
T01 = T1
T16 = T2*T3*T4*T5*T6
# pprint(simplify(T16.col(3)))
# print(latex(simplify(T16.col(3))))
# print('----')

# solveTheta1(T1, T16)
# solveTheta3(T01, T16)

# T03 = T1*T2*T3
# T36 = T4*T5*T6
# solveTheta2(T03,T36)

# T03 = T1*T2*T3
# T36 = T4*T5*T6
# solveTheta4(T03,T36)

# T04 = T1*T2*T3*T4
# T46 = T5*T6
# solveTheta5(T04,T46)

T05 = T1*T2*T3*T4*T5
T56 = T6
solveTheta6(T05,T56)


# T16a = T1.inv()*(T1*T2*T3*T4*T5*T6)
# T16b = T2*T3*T4*T5*T6
# cella =(T16a.col(3).row(1))
# cellb =(T16b.col(3).row(1))


# pprint(simplify(cellb))
