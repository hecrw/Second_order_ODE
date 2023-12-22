import cmath
import re
from sympy import symbols, integrate, exp, sin, cos, Mul, Matrix, diff, simplify, sympify, trigsimp
def assemble_case1(m1, m2):
    x, C1, C2 = symbols('x C1 C2')
    ycf = 0

    if m1 == 1.0:
        ycf += C1 * exp(x)
    else:
        ycf += C1 * exp(m1 * x)

    if m2 == 1.0:
        ycf += C2 * exp(x)
    else:
        ycf += C2 * exp(m2 * x)
    return ycf

def assemble_case2(m):
    x, C1, C2 = symbols('x C1 C2')
    ycf = 0

    if m == 1.0:
        ycf = exp(x) * C1 + exp(x) * C2 * x
    else:
        ycf =  exp(m*x) * C1 + exp(m*x) * C2 * x

    return ycf

def assemble_case3(m1):
    x, C1, C2 = symbols('x C1 C2')
    ycf = 0

    if m1.real == 0.0 and m1.imag == 1.0:
        ycf = C1 * sin(x) + C2 * cos(x)
    elif m1.real == 0.0:
        ycf = C1 * sin(m1.imag * x) + C2 * cos(m1.imag * x)
    elif m1.imag == 1.0:
        ycf =  exp(m1.real * x) * C1 * sin(x) + exp(m1.real * x) * C2 * cos(x)
    else:
        ycf =  exp(m1.real * x) * C1 * sin(m1.imag * x) + exp(m1.real * x) * C2 * cos(m1.imag * x)

    return ycf

def solve(ODE):
    
    if "y''" not in ODE: #not a second order ODE
        print("Enter a valid 2nd order ODE")
        return 0
    coeffecients = [1, 0, 0]
    fx = ""
    for i in range(len(ODE)): #this for loop in simple terms just collects a, b and c to solve the quadratic formula
        
        if ODE[i] == "y''": #skipping y''
            continue
        
        if ODE[i] == "-" or ODE[i] == "+":
            sign = ODE[i]
            continue
        
        if ODE[i] == "=":
            fx = ODE[i + 1:]
        if "y'" in ODE[i]:
            match = re.match(r"([-+]?\d*\.\d+|\d+)", ODE[i]) #extracting the coeffecient from the operand
            if match:
                coeffecients[1] = int(sign + str(match.group()))
            else:
                coeffecients[1] = int(sign + "1") #if no coeffecient exists then b will equal to 1
            continue
                
        elif "y" in ODE[i]:
            match = re.match(r"([-+]?\d*\.\d+|\d+)", ODE[i]) #extracting the coeffecient from the operand
            if match:
                coeffecients[2] = int(sign + str(match.group()))
            else:
                coeffecients[2] = int(sign + "1") #if no coeffecient exists then c will equal to 1
            continue
    
    a, b, c = coeffecients[0], coeffecients[1], coeffecients[2]
    root = b**2 - 4*a*c
    
    if root > 0: #first case 
        m1 = (-b + cmath.sqrt(root)) / (2 * a)
        m2 = (-b - cmath.sqrt(root)) / (2 * a)
        m1, m2 = m1.real, m2.real
        ycf = assemble_case1(m1, m2)
           
    elif root == 0: #second case
        m = (-b) / (2*a)
        m = m.real
        ycf = assemble_case2(m)
        
    elif root < 0:
        m1 = (-b + cmath.sqrt(root)) / (2 * a)
        ycf = assemble_case3(m1)
        
    if fx[0] == "0": #constant coeffecients
        return ycf
    
    else: #variation of parameters
        x, C1, C2 = symbols('x C1 C2')

        # Variation of parameters
        dummy = ycf.as_ordered_terms()
        functions = []
        for term in dummy:
            coefficient, function = term.as_coeff_mul()
            functions.append(function)
        y1 = [term for term in functions[0] if term != C1]
        y2 = [term for term in functions[1] if term != C2]

        y1 = Mul(*y1)
        y2 = Mul(*y2)
        
        mat = Matrix([[y1, diff(y1, x)], [y2, diff(y2, x)]])
        wronskian = simplify(mat.det())
        expression = "".join(fx)
        fx = sympify(expression)
        v1 = -integrate((y2 * fx) / wronskian)
        v2 = integrate((y1 * fx) / wronskian)
        yp = v1 * y1 + v2 * y2
        ygs = ycf + yp
        ygs = simplify(ygs)
        ygs = trigsimp(ygs)
        return ygs
        

ode = input("Enter your 2nd order ODE\n#Make sure each operand is seperated by a space\n")

solution = solve(ode.split())

print(solution)