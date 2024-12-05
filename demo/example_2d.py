import numpy as np



import polyquad

if __name__=='__main__':
    #pentagon
    verts = np.array(((1,-1), (-1,0), (0,3), (2,3), (3,0)))
    face = np.arange(5)

    k=4
    points, weights, r = polyquad.get_quadrature_2d(k, verts, face, get_residual = True)
    
