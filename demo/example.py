"""
In this example we show how to generate quadrature over an polyhedron
"""
from time import perf_counter

from scipy.io import loadmat

import polyquad

if __name__=='__main__':
    case = 'icosphere'
    geo = loadmat(case + '.mat')
    # get the coordinates of vertices
    verts = geo['verts']
    # get the connectivity table
    faces = geo['faces']

    integrand_order = 20
    points, weights, residual = polyquad.get_quadrature_3d(integrand_order, verts, faces, get_residual = True)
    #NOTE: the first time you run the script it may take a bit of time, don't panic!
    #There are two reasons for that:
    #     - some functions are wrapped with numba, so the first time you call them some compiling is happening. This is done only once per python session
    #     - a QR decomposition is performed at some point. Luckily, this decomposition only depends on the polynomial order and can be reused for different shapes.

    t1 = perf_counter()
    _,_ = polyquad.get_quadrature_3d(19, verts, faces)
    t2 = perf_counter()
    _,_ = polyquad.get_quadrature_3d(19, verts, faces)
    t3 = perf_counter()
    print(f"Elapsed time to generate quadrature of order 17 for the {case}: with computation of QR:\n\t{t2-t1}s")
    print(f"Elapsed time to generate quadrature of order 17 for the {case}: reusing the previously computer QR:\n\t{t3-t2}s")
