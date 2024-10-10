import numpy as np
import numba

def trios(order: int) -> (np.ndarray):
    """combination of powers of monomials for a polynomial of a given order with complete basis

    Parameters
    ----------
    order : int
        polynomial order

    Returns
    -------
    xPow np.ndarray
        powers for x
    yPow np.ndarray
        powers for y
    zPow np.ndarray
        powers for z
    """
    jj = np.arange(0,order+1, dtype = 'int32')
    v1,v2,v3 = np.meshgrid(jj,jj,jj)
    mask = v1+v2+v3 < order+1
    yPow, xPow, zPow = v1[mask], v2[mask], v3[mask]
    return xPow.reshape((xPow.size,1)), yPow.reshape((yPow.size,1)),zPow.reshape((zPow.size,1))

def duos(order: int) -> (np.ndarray):
    """combination of powers of monomials for a polynomial of a given order with complete basis

    Parameters
    ----------
    order : int
        polynomial order

    Returns
    -------
    xPow np.ndarray
        powers for x
    yPow np.ndarray
        powers for y
    """
    jj = np.arange(0,order+1, dtype = 'int32')
    v1,v2 = np.meshgrid(jj,jj)
    mask = v1+v2 < order+1
    xPow, yPow = v1[mask], v2[mask]
    return yPow.reshape((yPow.size,1)), xPow.reshape((xPow.size,1))


@numba.jit(nopython = True, cache = True, nogil = True)
def getFaceNormal(face: np.ndarray,
               verts: np.ndarray):
    #https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
    numEdges = face.shape[0]
    normal = np.array((0.,0.,0.))
    for ii in range(numEdges):
        currentVert = verts[face[ii],:]
        nextVert = verts[face[(ii+1) % numEdges],:]
        normal[0] += (currentVert[1] - nextVert[1]) * (currentVert[2] + nextVert[2])
        normal[1] += (currentVert[2] - nextVert[2]) * (currentVert[0] + nextVert[0])
        normal[2] += (currentVert[0] - nextVert[0]) * (currentVert[1] + nextVert[1])
    normal = normal / np.sqrt(np.sum(normal**2))
    return normal
