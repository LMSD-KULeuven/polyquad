import numpy as np

import polyquad

if __name__=='__main__':
    #pentagon
    verts = np.array(((1,-1), (-1,0), (0,3), (2,3), (3,0)))
    face = np.arange(5)

    k=4
    points, weights, r = polyquad.get_quadrature_2d(k, verts, face, get_residual = True)
    

    #polygon 7.d fron chin and lassere
    x = [1,3,3,1]
    y = [-.5,-.5,1.5,1.5]
    pts = np.vstack((x,y)).T
    face = np.arange(5)
    p,w = polyquad.get_quadrature_2d(20,pts,face)

    # x= p[:,0]
    # y= p[:,1]
    # I = sum((x**2 + x*y + y**2)*w)

    p,j = polyquad.map_to_local_bb_2d(pts)
    I = polyquad.integrateMonomialsPolygon(4,np.array([0,1,2,3]),np.pad(p,((0,0),(0,1))))
    p,w = polyquad.get_quadrature_2d(4 , p, np.arange(4))
