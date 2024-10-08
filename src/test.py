import time

import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt

import polyquad


if __name__=='__main__':
    geometry_file = '../demo/bunny.mat'
    geo = loadmat(geometry_file)
    vertices = geo['verts']
    vertices[:,0] = vertices[:,0]*2
    vertices +=3
    faces = geo['faces']

    orders = np.arange(1,25)
    resList = np.zeros(len(orders))
    resListInv = np.zeros(len(orders))
    for ii,order in enumerate(orders):
        print(f'doing order {order}')
        p, w, r = polyquad.get_quadrature_3d(order, vertices, faces, get_residual = True)
        resList[ii] = r

    plt.semilogy(orders, resList, label='solve(a,b)')
    plt.legend()
    plt.grid()
    plt.show()



    
