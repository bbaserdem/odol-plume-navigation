#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MYSPH : SPH simulation for odorant plume simulation
    Copyright (C) 2019  Batuhan Başerdem

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import sys

import numpy as n

def init_particle_box(num=5, size=.1, dim=2, pack=2):
    """ Generate particles equally spaced in dim dimensional box """
    # Generator for introducing new dimensions
    pgen = n.array([
        n.linspace(
            (size/2)*(pack+1)/(num*pack+1),
            (size/2)*((2*num-1)*pack+1)/(num*pack+1),
            num)
        ]).T
    # Generator for capping new dimensions with boundaries
    bgen = n.array([n.linspace(0, size, pack*num+2)]).T
    # Points and boundary points
    pts = n.empty((1, 0))
    bnd = n.empty((0, 0))
    # Storing the caps to append to the boundary
    cgen = n.array([n.linspace(0, size, 2)]).T
    cap = n.empty((1, 0))
    for _ in range(dim):
        # Multiply existing points along the new axis
        pts = n.append(n.kron(n.ones((num, 1)), pts),
                       n.kron(pgen, n.ones((pts.shape[0], 1))),
                       axis=1)
        # Multiply boundary accross the new dimension, including the edge
        bnd = n.append(n.kron(n.ones((pack*num+2, 1)), bnd),
                       n.kron(bgen, n.ones((bnd.shape[0], 1))),
                       axis=1)
        # Append the cap at the maxima of new dimension
        tip = n.append(n.kron(n.ones((2,1)), cap),
                       n.kron(cgen, n.ones((cap.shape[0], 1))),
                       axis=1)
        bnd = n.append(bnd,
                       tip,
                       axis=0)
        # Create the new cap
        cap = n.append(n.kron(n.ones((pack*num, 1)), cap),
                       n.kron(bgen[1:-1,:], n.ones((cap.shape[0], 1))),
                       axis=1)
    print(pgen)
    print(pts)
    print(bnd)

init_particle_box(num=2, size=1, dim=3)
