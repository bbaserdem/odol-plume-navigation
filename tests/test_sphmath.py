#!/usr/bin/python
"""
    Test methods for the math module of SPH
"""


# Builtin modules

# Third party modules
import numpy as n
import scipy.sparse as s

# Self modules
from odorseeker.sph import math as m

bord1d = n.array([[0], [1]])
bord1d = n.array([
    [0, ],[1]])

ex1d = {
        'particles': n.array([[.5]]),
        :q
