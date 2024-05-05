__author__ = 'mkv-aql'
from sympy import Symbol, cos, factor, symbols

x = Symbol('x')
e = 1/cos(x)
print(e.series(x, 0, 10))

print(x+x+1)

y = Symbol('y')

#Factorization
square = x**2 - 2*x*y + y**2
res = factor(square)
print(res)

#Simple declaration
x,y,z = symbols('x,y,z')
print(x)
print(y)
print((z*z)/y)

print((z*z)/z)
print(z*2 + x*0)


#ODE solving
from sympy import Function, dsolve, Eq, Derivative
f = Function('f')
result  = dsolve(Derivative(f(x), x) - f(x), f(x))
print(result)

result_1 = Derivative(f(x), x)
print(result_1)
print(dsolve(result_1))


#Integration
from sympy import integrate, init_printing
k, j = symbols('k,j')
d = integrate(k, (j, a, b))
print(d)
d = integrate(k, (j, -1, 1))
print(d)
init_printing(use_unicode=False, wrap_line=False)
d = integrate(k**2 + k + 1, k)
print(d)