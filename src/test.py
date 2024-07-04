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

    order = 18
    p, w, r = polyquad.get_quadrature(order, vertices, faces, get_residual = True)
    p, w = polyquad.get_quadrature(order, vertices, faces, get_residual = False)

    ntests = 10
    t1 = time.perf_counter()
    for _ in range(ntests):
        p, w = polyquad.get_quadrature(order, vertices, faces, get_residual = False)
    t2 = time.perf_counter()
    print(f'elasped time {(t2-t1)/ntests}s')
    
    v2,j = polyquad.map_to_local_bb(vertices)
    I = polyquad.integrateMonomialsPolyhedron(10,faces,v2)
