"""A python wrapper around fast C++ routines for numemrically integrating ODEs."""

from ._pyodepack import solveIVP
from . import methods


__all__ = ['solveIVP', 'methods']
