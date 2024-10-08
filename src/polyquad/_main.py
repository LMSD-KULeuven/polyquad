import numpy as np

from ._mapping import map_to_local_bb
from ._antonietti import integrateMonomialsPolyhedron
from ._moment_matching import moment_matching

def get_quadrature(order:int,
                   vertices: np.ndarray,
                   faces: np.ndarray,
                   get_residual = False)-> (np.ndarray, np.ndarray):
    """generates quadrature points and weights for a specified order

    Parameters
    ----------
    order : int
        total polynomial order of the integrand
    vertices : np.ndarray
        coordinates of vertices, shape = (numver of vertices, 3)
    faces : np.ndarray
        declaration of planar faces with ordered vertices indices. If
        all faces have the save number of vertices then "faces" can be
        a np.ndarray of shape (number of faces, number of vertices per
        face)
    get_residual : bool
        if true returns the residual as defined in Eq.17 of the paper.
        It is the direct measure of the L2-norm error on the
        integration of each monomial.

    Returns
    -------
    (np.ndarray, np.ndarray)
        points, weights

    """
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
