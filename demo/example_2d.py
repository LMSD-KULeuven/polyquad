"""
In this example we show how to generate quadrature over a polygon
"""
from time import perf_counter

import polyquad
import numpy as np


if __name__=='__main__':
    verts = np.array(((1,-1), (-1,0), (0,3), (2,3), (3,0)))
    face = np.array((0,1,2,3,4))

    k = 4

    #the polygon doesn't lie in [-1,1]^2 so we need to ask for a mapping
    points, weights = polyquad.get_quadrature_2d(k, verts, face, mapping = True)
    
    print(sum(weights))
