from __future__ import print_function, division
from sympy import *
from sympy.utilities import lambdify
from scipy.integrate import ode as scipy_ode
import numpy as n


def extract_variables(expression):
    """Finds the (symbolic) variables in an expression"""
    
    if 'matrices' in str(type(expression)):
        args = list(expression)
    elif len(expression.args) == 0:
        return {expression}
    else:
        args = expression.args

    variables = set()

    for var in args:
        variables.update(extract_variables(var))
        
    return set([var for var in variables if 'number' not in str(type(var))])





