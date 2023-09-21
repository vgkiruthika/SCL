#ps-3

#bisection method
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import sys


def func(eq, val):
  return eq.subs(x,val)

def bisection(eq,a,b):
  if (func(eq,a)*func(eq,b)) >=0:
    print("Invalid a and b")

  c=a
  while((b-a)>=0.01):
    c=(a+b)/2

    if(func(eq,c)==0):
      break

    if(func(eq,c)*func(eq,a)<0):
      b=c

    elif(func(eq,c)*func(eq,b)<0):
      a=c

    print("{0:.4f}\t\t{1:.4f}\t\t{2:.4f}\t\t{3:.4f}\t\t{4:.4f}\t\t{5:.4f}".format(a,b,c,func(eq,c),func(eq,c)*func(eq,a),func(eq,c)*func(eq,b)))


  return c

x = symbols('x')
inp = input("Enter the function:")
equation = Eq(eval(inp), 0)
print(equation)
a = int(input("Enter a:"))
b = int(input("Enter b:"))
print(equation)
lhs = equation.lhs

# Plotting the function
p=plot(equation.lhs,(x,a,b))

#find f(a) and f(b)
f_a=func(lhs,a)
f_b=func(lhs,b)
if(f_a==0):
  print("Root: ",a)
  sys.exit()
elif(f_b==0):
  print("Root: ",b)
  sys.exit()

fa_fb=f_a*f_b
if(fa_fb>0):
  print("Inappropriate interval!")
  sys.exit()

print("a\t\tb\t\tmid\t\tf(mid)\t\tf(a)*f(mid)\tf(b)*f(mid)")

root=bisection(lhs,a,b)
print("\n\nROOT: ", "%.4f"%root)



#Regula falsi method
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import sys


def func(eq, val):
  return eq.subs(x,val)

def regula_falsi(eq,a,b):
  h=((abs(func(eq,a)))*(b-a))/(abs(func(eq,a))+abs(func(eq,b)))
  c=a+h
  print("{0:.4f}\t\t{1:.4f}\t\t{2:.4f}\t\t{3:.4f}\t\t{4:.4f}\t\t{5:.4f}\t\t{6:.4f}".format(a,b,func(eq,a),func(eq,b),h,c,func(eq,c)))
  while True:

     if func(eq,c)==0:
      break
     temp=c
     if(func(eq,c)<0):
       a=c
     else:
       b=c
     h=((abs(func(eq,a)))*(b-a))/(abs(func(eq,a))+abs(func(eq,b)))
     c=a+h
     print("{0:.4f}\t\t{1:.4f}\t\t{2:.4f}\t\t{3:.4f}\t\t{4:.4f}\t\t{5:.4f}\t\t{6:.4f}".format(a,b,func(eq,a),func(eq,b),h,c,func(eq,c)))
     if abs(temp-c)<=0.00001 :
       break
  return c

x = symbols('x')
inp = input("Enter the function:")
equation = Eq(eval(inp), 0)
print(equation)
a = int(input("Enter a:"))
b = int(input("Enter b:"))
print(equation)
lhs = equation.lhs

# Plotting the function
p=plot(equation.lhs,(x,a,b))

#find f(a) and f(b)
f_a=func(lhs,a)
f_b=func(lhs,b)
if(f_a==0):
  print("Root: ",a)
  sys.exit()
elif(f_b==0):
  print("Root: ",b)
  sys.exit()

fa_fb=f_a*f_b
if(fa_fb>0):
  print("Inappropriate interval!")
  sys.exit()
print("a\t\tb\t\tf(a)\t\tf(b)\t\th\t\tx\t\tf(x)")
root=regula_falsi(lhs,a,b)
print("\n\nROOT: ", "%.4f"%root)



#newton raphson method

from sympy import *
import matplotlib.pyplot as plt

x = symbols('x')
inp = input("Enter the function:")
equation = Eq(eval(inp), 0)
print(equation)


a=float(input("Enter starting value of interval:"))
b=float(input("Enter ending value of interval:"))
fun=equation.lhs
deriv=diff(fun,x)
p=plot(equation.lhs,(x,a,b))
x_val=a
print("x_n\t\t\tf(xn)\t\t\tf'(xn)\t\t\tx_n+1")

while(True):
  temp=x_val
  x_val=temp-(fun.subs(x,temp)/deriv.subs(x,temp))
  print("{0:.6f}\t\t{1:.6f}\t\t{2:.6f}\t\t{3:.8f}".format(temp,fun.subs(x,temp),deriv.subs(x,temp),x_val))
  if(abs(temp-x_val)<0.0000001):
    break

print("\n\nSolution: {0:.4f}".format(temp))


#fixed point iteration method
from sympy import *
import matplotlib.pyplot as plt

import sys
x = symbols('x')
inp = input("Enter the function f(x):")
f_eq = Eq(eval(inp), 0)
print(f_eq)
f=f_eq.lhs
inp=input("Enter the function phi(x)")
phi_eq=Eq(eval(inp),0)
phi=phi_eq.lhs

a=float(input("Enter starting value of interval:"))
b=float(input("Enter ending value of interval:"))
p=plot(f,(x,a,b))
deriv_phi=diff(phi)
d_a=abs(deriv_phi.subs(x,a))
d_b=abs(deriv_phi.subs(x,b))
maxi=max(d_a,d_b)

if(maxi<1):
  print("phi(x) gives us a convergent sequence of iteration.")
else:
  print("Not convergent!")
  sys.exit()

x_val=a
print("x_n\t\t\tphi(x_n)")
while(True):
  temp=x_val
  x_val=phi.subs(x,temp)
  print("{0:.6f}\t\t{1:.6f}".format(temp,x_val))
  if(abs(temp-x_val)<0.00001):
    break

print("\n\nSolution: {0:.4f}".format(temp))

#ps-4

#1) Linear interpolation

from sympy import *
import matplotlib.pyplot as plt

x,y=symbols('x y')
x1=float(input("Enter x1:"))
y1=float(input("Enter y1:"))
x2=float(input("Enter x2:"))
y2=float(input("Enter y2:"))

y=(y1+((y2-y1)/(x2-x1))*(x-x1))
print(y)

inp=float(input("Enter the point of interpolation:"))
res=y.subs(x,inp)
print("Res: ",res)

x_vals=[x1,x2,inp]
y_vals=[]
for val in x_vals:
  y_vals.append(y.subs(x,val))

fig, ax = plt.subplots()
ax.plot(x_vals, y_vals)

ax.scatter(x_vals, y_vals, color='blue')


#2)Quadratic interpolation

from sympy import *
import matplotlib.pyplot as plt
import numpy as np

x,y=symbols('x y')
x0=float(input("Enter x0:"))
y0=float(input("Enter y0:"))
x1=float(input("Enter x1:"))
y1=float(input("Enter y1:"))
x2=float(input("Enter x2:"))
y2=float(input("Enter y2:"))


b0=y0
b1=(y1-y0)/(x1-x0)
b2=(((y2-y1)/(x2-x1))-b1)/(x2-x0)
print(b0)
print(b1)
print(b2)

y=b0+b1*(x-x0)+(b2*(x-x0)*(x-x1))


print(y)

inp=float(input("Enter the point of interpolation:"))
res=y.subs(x,inp)
print("Res: ",res)
x_pts=[x0,x1,x2]
y_pts=[]

for pt in x_pts:
  y_pts.append(y.subs(x,pt))

x_vals=np.linspace(x0,x2,100)
y_vals=[]
for val in x_vals:
  y_vals.append(y.subs(x,val))

fig, ax = plt.subplots()
ax.plot(x_vals, y_vals)

ax.scatter(x_pts, y_pts, color='blue')


#3)Quadratic interpolation

from sympy import *
import matplotlib.pyplot as plt
import numpy as np

x,y=symbols('x y')
x0=0
y0=0
x1=1
y1=1
x2=2
y2=20


b0=y0
b1=(y1-y0)/(x1-x0)
b2=(((y2-y1)/(x2-x1))-b1)/(x2-x0)
print(b0)
print(b1)
print(b2)

y=b0+b1*(x-x0)+(b2*(x-x0)*(x-x1))


print(y)

x_pts=[x0,x1,x2]
y_pts=[]

for pt in x_pts:
  y_pts.append(y.subs(x,pt))

x_vals=np.linspace(x0,x2,100)
y_vals=[]
for val in x_vals:
  y_vals.append(y.subs(x,val))

fig, ax = plt.subplots()
ax.plot(x_vals, y_vals)

ax.scatter(x_pts, y_pts, color='blue')

#4)Newton's divided difference
from sympy import *

def divided_diff(x_val, y_val, table, n):
    for i in range(1, n):
        for j in range(n - i):
            table[j][i] = (table[j+1][i-1] - table[j][i-1]) / (x_val[i+j] - x_val[j])
    return table

n = int(input("Enter the number of points: "))
x, y = symbols('x y')
x_val = []
y_val = []

for i in range(n):
    inpx = float(input("Enter x{0}: ".format(i)))
    inpy = float(input("Enter y{0}: ".format(i)))
    x_val.append(inpx)
    y_val.append(inpy)

table = [[0 for i in range(n)] for i in range(n)]

for i in range(n):
    table[i][0] = y_val[i]

table = divided_diff(x_val, y_val, table, n)

print("\n\nDivide and diference table:")
for row in table:
    print(row)

interpolation_term = 0
for i in range(n):
    term = table[0][i]
    for j in range(i):
        term *= (x - x_val[j])
    interpolation_term += term

result = simplify(interpolation_term)
print("\n\nInterpolating polynomial:\n\n", result)

inp = float(input("Enter the interpolated value: "))
res = result.subs(x, inp)
print("\n\nInterpolated result:", res)

x_pts=np.linspace(-5,5,100)
y_pts=[]
for val in x_pts:
  y_pts.append(result.subs(x,val))

fig, ax = plt.subplots()
ax.plot(x_pts, y_pts)

ax.scatter(x_val, y_val, color='blue')


#5)Newton's divided difference
from sympy import *

def divided_diff(x_val, y_val, table, n):
    for i in range(1, n):
        for j in range(n - i):
            table[j][i] = (table[j+1][i-1] - table[j][i-1]) / (x_val[i+j] - x_val[j])
    return table

n = int(input("Enter the number of people: "))
x, y = symbols('x y')
x_val = []
y_val = []

for i in range(n):
    inpx = float(input("Enter age{0}: ".format(i)))
    inpy = float(input("Enter premium{0}: ".format(i)))
    x_val.append(inpx)
    y_val.append(inpy)

table = [[0 for i in range(n)] for i in range(n)]

for i in range(n):
    table[i][0] = y_val[i]

table = divided_diff(x_val, y_val, table, n)

print("\n\nDivide and difference table:")
for row in table:
    print(row)

interpolation_term = 0
for i in range(n):
    term = table[0][i]
    for j in range(i):
        term *= (x - x_val[j])
    interpolation_term += term

result = simplify(interpolation_term)
print("\n\nInterpolating polynomial:\n\n", result)

inp = float(input("Enter the interpolated value: "))
res = result.subs(x, inp)
print("\n\nInterpolated result:", res)

x_pts=np.linspace(x_val[0]-10,x_val[n-1]+10,100)
y_pts=[]
for val in x_pts:
  y_pts.append(result.subs(x,val))

fig, ax = plt.subplots()
ax.plot(x_pts, y_pts)
x_val.append(inp)
y_val.append(res)
ax.scatter(x_val, y_val, color='blue')


#ps-5

#Lagrange’s interpolation method
import matplotlib.pyplot as plt
import numpy as np

print("Enter x values:")
x=list(map(float,input().split()))
print("Enter f(X) values:")
f_x=list(map(float,input().split()))
x1=float(input("Enter x value to find f(x):"))

value =0

for i in range(len(x)):
    num=1
    den=1
    for k in x:
      if k!=x[i]:
        num*=(x1-k)
        den*=(x[i]-k)
    value+=((num/den)*f_x[i])

print("f({0}) value is {1:.4f}".format(x1,value))

x_values = np.linspace(min(x), max(x), 100)
y_values = []

for x_val in x_values:
    y = 0
    for i in range(len(x)):
        num = 1
        den = 1
        for k in range(len(x)):
            if k != i:
                num *= (x_val - x[k])
                den *= (x[i] - x[k])
        y += ((num / den) * f_x[i])
    y_values.append(y)

plt.plot(x_values, y_values, label='Interpolated Polynomial', linestyle='-', color='b')
plt.scatter(x, f_x, c='blue', label='Data Points')  # Scatter plot for data points
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title("Lagrange's Interpolation")
plt.legend()
plt.grid(True)
plt.show()

#Inverse lagrange’s interpolation method
print("Enter x values:")
x=list(map(float,input().split()))
print("Enter f(X) values:")
f_x=list(map(float,input().split()))
y=float(input("Enter f(x) value to find x:"))

value =0

for i in range(len(x)):
    num=1
    den=1
    for k in f_x:
      if k!=f_x[i]:
        num*=(y-k)
        den*=(f_x[i]-k)
    value+=((num/den)*x[i])

print("f-1({0}) value is {1:.4f}".format(y,value))

plt.scatter(f_x, x, c='blue', label='Data Points')  # Scatter plot for data points
plt.axhline(y=value, color='blue', linestyle='--', label=f'x = {value}')  # Line for the calculated x
plt.xlabel('f(x)')
plt.ylabel('x')
plt.title("Inverse Lagrange's Interpolation")
plt.legend()
plt.grid(True)
plt.show()


#Cubic spline interpolation

import sympy as sp
from sympy import *
from math import *
import matplotlib.pyplot as plt
import numpy as np

n = int(input("Enter the number of data points:"))
x_pts=list(map(int,input("Enter x vals:").split()))
y_pts=list(map(int,input("Enter y vals:").split()))

m = []
for i in range(n):
    m.append(Symbol(f'm_{i}'))
print(m)

equations=[]

for i in range(1,n-1):=
  print("in")
  eq=Eq((m[i-1]+(4*m[i])+m[i+1]),6*(y_pts[i-1]-(2*y_pts[i])+y_pts[i+1]))
  equations.append(eq)

substitutions = {m[0]: 0, m[n-1]: 0}
equations_with_subs = [eq.subs(substitutions) for eq in equations]

m_final=[]
rem_m=[]
for i in range(0,n):
  if i==0 or i==n-1:
    m_final.append(0)
  else:
    m_final.append(m[i])
    rem_m.append(m[i])


for eq in equations_with_subs:
  print(eq)

solution = sp.solve(equations_with_subs,rem_m)
print(solution)

x=Symbol('x')
f_x=[]
for i in range(0,n-1):
  print(i)
  y=(1/6)*((x_pts[i+1]-x)**3)*m_final[i] + (1/6)*((x-x_pts[i])**3)*m_final[i+1] + (x_pts[i+1]-x)*(y_pts[i]-(1/6)*m_final[i]) + (x-x_pts[i])*(y_pts[i+1]-(1/6)*m_final[i+1])
  f_x.append(y)

res_f=[]

for f in f_x:
    for i in range(1, n - 1):
        f = f.subs(sp.Symbol(f'm_{i}'), solution[sp.Symbol(f'm_{i}')])
    res_f.append(simplify(f))
    print(simplify(f))


inp=float(input("Enter the interpolated value:"))
c_inp=ceil(inp)
res=res_f[int(c_inp)-1].subs(x,inp)
print("Result:",res)

x_vals_all = []
y_vals_all = []

for i in range(n-1):
  x_vals=[]
  y_vals=[]
  x1=x_pts[i+0]
  x2=x_pts[i+1]
  x_vals.extend(np.linspace(x1,x2,50))
  for val in x_vals:
    y_vals.append(res_f[i].subs(x,val))
  x_vals_all.extend(x_vals)
  y_vals_all.extend(y_vals)

fig, ax = plt.subplots()

x_pts.append(inp)
y_pts.append(res)
ax.scatter(x_pts, y_pts, color='blue', label='Data Points')
ax.plot(x_vals_all, y_vals_all, label='Interpolation')

plt.legend()
plt.show()


#input:
#x:1 2 3 4
#y:1 2 5 11
#inp:1.5
#res: 1.375

#NEWTON'S FORWARD INTERPOLATION

import matplotlib.pyplot as plt
def create_table(x,y,n,table):
  for col in range(1,n):
    for row in range(n-col):
      table[row][col]=table[row+1][col-1]-table[row][col-1]

  return table

def calculate(u,n):
  temp = u;
  for i in range(1, n):
      temp = temp * (u - i);
  return temp;

def fact(n):
    f = 1;
    for i in range(2, n + 1):
        f *= i;
    return f;

x=list(map(float,input("Enter x values:").split()))
y=list(map(float,input("Enter y values:").split()))

n=len(x)

table=[[0.0 for i in range(n)]for i in range(n)]
for i in range(n):
    table[i][0]=y[i]
#creating backward difference
table=create_table(x,y,n,table)
for row in table:
  print(row)


inp=float(input("Enter the value to be interpolated at:"))

sum=table[0][0]
u=(inp-x[0])/(x[1]-x[0])
polynomial = f"{y[0]:.2f}"

for i in range(1,n):
  sum = sum + (((calculate(u, i)) * (table[0][i])) / fact(i));
  polynomial += f" + {table[0][i]:.2f} * u"


print("Result: ",sum)
print("Interpolated Polynomial:", polynomial)

x_values = [x[0] + i * 0.01 for i in range(int((x[n-1] - x[0]) / 0.01) + 1)]
y_values = []


for val in x_values:
    u = (val - x[0]) / (x[1] - x[0])
    interpolated_y = table[0][0]
    for i in range(1, n):
        interpolated_y += (calculate(u, i) * table[0][i]) / fact(i)
    y_values.append(interpolated_y)

x.append(inp)
y.append(sum)

plt.plot(x_values, y_values, label='Interpolated Polynomial', linestyle='-', color='b')
plt.scatter(x, y, c='blue', label='Data Points')
plt.ylabel('y')
plt.title('Newton\'s Forward Interpolation')
plt.legend()
plt.grid(True)

#40 50 60 70 80
#31 73 124 159 190