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

import numpy as n
import scipy.sparse as s

# Error handling
class Error(Exception):
    """ Base class for errors """
    pass

class VectorMatrixInvalidInputError(Error):
    """ Raise this error when class encounters dimensional mismatch """
    pass

class VectorMatrixMismatchDimError(Error):
    """ Raise this error when class encounters dimensional mismatch """
    pass

class VectorMatrixNonSparseError(Error):
    """ Raise this error when class deals with non-sparse values """
    pass

# Class to create vector matrices
class VectorMatrix():
    """
    Class to store matrix of vectors
    """
    def __init__(self, *args):
        """
        Constructor
        Take list of sparse matrices, and merge them in one
        """
        self.data = args
        self.dim = len(args)
        self.shape = data[0].shape
        self.dtype = data[0].dtype
        # Checks
        for i in range(self.dim):
            # Shape check
            if not cmp(self.data[i].shape, self.shape):
                raise VectorMatrixInvalidInputError
            # Type check
            elif self.data[i].dtype != self.dtype:
                raise VectorMatrixInvalidInputError
            # Sparseness check
            elif not isinstance(self.data[i], s.spmatrix):
                raise VectorMatrixInvalidInputError

    def transpose(self):
        """ Overload the scipy transpose operator, apply to each spmatrix """
        out = list()
        for i in range(self.dim):
            out.append(self.data[i].T)
        return VectorMatrix(*out)

    def __matmul__(self, mat):
        """ Overloading of @ operator for these matrices """
        # SELF MULTIPLY
        if isinstance(mat, SparseVector):
            # Can only do it if both have the same dimensions
            if self.dim == mat.dim:
                if self.shape[1] == mat.shape[0]:
                    out = list()
                    for i in range(self.dim):
                        out.append(self.data[i].dot(mat.data[i]))
                    return VectorMatrix(*out)
                else:
                    raise VectorMatrixMismatchError
            else:
                raise VectorMatrixMismatchDimError
        # NUMPY ARRAY or SPARSEARRAY MULTIPLY
        elif isinstance(mat, n.ndarray) or isinstance(mat, s.spmatrix):
            # IF IS VECTOR
            if len(mat.shape) == 1:
                # Default: Contract by particle index
                if  self.shape[1] == mat.shape[0]:
                    out = n.empty((self.shape[0], self.dim))
                    for i in range(self.dim):
                        out[:,i] = self.data[i].dot(mat)
                    return out
                # Allow for contracting by spatial index
                elif mat.shape[0] == self.dim:
                    out = s.csc_matrix(self.shape, dtype=self.dtype)
                    for i in range(self.dim):
                        out += mat[i]*self.data[i]
                    return out
                else:
                    raise VectorMatrixMismatchDimError
            # IF IS MATRIX
            elif len(mat.shape) == 2:
                # Default: Matrices are multiplied on particle index
                if self.shape[1] == mat.shape[0]:
                    out = list()
                    for i in range(self.dim):
                        out.append(s.csc_matrix(self.data[i].dot(mat)))
                    return VectorMatrix(*out)
                # Also allow for contracting with spatial index
                elif self.dim == mat.shape[0]:
                    out = list()
                    for i in range(self.dim):
                        out.append(
                            s.csc_matrix(
                                self.shape,
                                dtype=self.dtype,
                                nnz=self.data[i].nnz))
                    for i in range(mat.shape[0]):
                        for j in range(mat.shape[1]):
                            out[j] += mat.data[i] * mat[i,j]
                    return VectorMatrix(*out)
                else:
                    raise VectorMatrixInvalidInputError
            # DID NOT CONFIGURE CONTRACTION WITH MULTIDIM
            else:
                return NotImplemented
        # IF IS SPARSE
        else:
            return NotImplemented

    def __mul__(self, mat):
        """ Overloading of * operator for these matrices """
        # SELF PRODUCT
        if isinstance(mat, SparseVector):
            # Can only do it if both have the same dimensions
            if self.dim == mat.dim:
                if self.shape == mat.shape:
                    out = list()
                    for i in range(self.dim):
                        out.append(self.data[i].multiply(mat.data[i]))
                    return VectorMatrix(*out)
                else:
                    raise VectorMatrixMismatchError
            else:
                raise VectorMatrixMismatchDimError
        # NUMPY ARRAY or SPARSEARRAY MULTIPLY
        elif isinstance(mat, n.ndarray) or isinstance(mat, s.spmatrix):
            # IF IS VECTOR
            if len(mat.shape) == 1:
                # Default: Multiply at particle index
                if self.shape[1] == mat.shape[0]
                    out = list()
                    for i in range(self.dim):
                        out.append(self.data[i].multiply(mat))
                    return VectorMatrix(*out)
                # Allow for multiplication with spatial too
                elif self.dim == mat.shape[0]:
                    out = list()
                    for i in range(self.dim):
                        out.append(mat[i]*self.data[i])
                    return VectorMatrix(*out)
                else:
                    raise VectorMatrixMismatchDimError
            # IF IS MATRIX
            elif len(mat.shape) == 2:
                # Default: Matrices are multiplied on particle index
                if self.shape == mat.shape:
                    out = list()
                    for i in range(self.dim):
                        out.append(self.data[i].multiply(mat))
                    return VectorMatrix(*out)
                else:
                    raise VectorMatrixInvalidInputError
            # DID NOT CONFIGURE MULTIDIM
            else:
                return NotImplemented
        # If semething else, try until fail
        else:
            out = list()
            for i in range(self.dim)A
                out.append(self.data[i]*mat)
            return VectorMatrix(*out)
    
    def __rmatmul__(self, mat):
        """ Reverse operation of matrix multiplication, abuse transposes """
        out = self.transpose().__matmul__(mat.transpose())
        return out.transpose()
    
    def __rmatmul__(self, mat):
        """ Reverse operation of matrix product, abuse transposes """
        out = self.transpose().__mul__(mat.transpose())
        return out.transpose()

def distance_all_by_all(cor, thr):
    """ Distance function to get pairwise distances under threshold """

if __name__ == "__main__":
    # Do a display of initialization if called explicitly
    
