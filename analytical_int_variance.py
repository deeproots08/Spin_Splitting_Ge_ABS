# Import sympy for symbolic computation
import numpy as np
import sympy as sp

# Redefine variables
x, y = sp.symbols('x y', real=True)

# Define the expressions for each case
expr1 = sp.exp(-(x - y)) * sp.sin(x * sp.pi / 9)**2 * sp.sin(y * sp.pi / 4)**2
expr2 = sp.exp(-(y - x)) * sp.sin(x * sp.pi / 9)**2 * sp.sin(y * sp.pi / 4)**2

# Set up and compute the integrals for each case
# Integral for case 1: x ranges from 0 to 9, y ranges from 0 to x
integral_case1 = sp.integrate(expr1, (y, 0, x), (x, 0, 9))
# Integral for case 2: x ranges from 0 to 4, y ranges from x to 4
integral_case2 = sp.integrate(expr2, (y, x, 4), (x, 0, 4))

# Total integral
total_integral = integral_case1 + integral_case2
total_integral.evalf()  # Evaluate to get a numerical approximation if necessary
