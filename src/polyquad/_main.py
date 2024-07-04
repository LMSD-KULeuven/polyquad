import numpy as np

from ._mapping import map_to_local_bb
from ._antonietti import integrateMonomialsPolyhedron
from ._moment_matching import moment_matching

def get_quadrature(order:int,
                   vertices: np.ndarray,
                   faces: np.ndarray,
                   get_residual = False):
    # perform matrix transfomation to fit in the local bounding box = [-1, 1]^3
    verts, jacobian = map_to_local_bb(vertices)
    # call antonietti's algorithm
    integrated_monomials = integrateMonomialsPolyhedron(order, faces, verts)
    # perform moment matching
    if get_residual:
        points, weights, resiadual = moment_matching(order, integrated_monomials, residual = True)
        weights = weights*jacobian
        return points, weights, resiadual
    else:
        points, weights = moment_matching(order, integrated_monomials, residual = False)
        weights = weights*jacobian
        return points, weights
