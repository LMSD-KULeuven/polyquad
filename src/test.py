import time

from scipy.io import loadmat

import polyquad

if __name__=='__main__':
    geometry_file = '../demo/icosphere.mat'
    geo = loadmat(geometry_file)
    vertices = geo['verts']
    vertices[:,0] = vertices[:,0]*2
    vertices +=3
    faces = geo['faces']

    order = 10

    p, w, r = polyquad.get_quadrature(order, vertices, faces, get_residual = True)
