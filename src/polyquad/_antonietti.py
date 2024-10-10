"""
Implementation of Paola Antonietti's algorithm 2
https://doi.org/10.1007/s10915-018-0802-y

Hats off to Thijs van Putten coding it in the first place.
"""

import numpy as np
import numba

from ._utils import getFaceNormal

@numba.jit(nopython = True, cache = True, nogil = True)
def getFaceCenter(verts: np.ndarray) -> np.ndarray:
    """computes the barycenter of vertices

    Could be done with a simple mean, but unavailable in numba

    Parameters
    ----------
    verts : np.ndarray
        coordinates of vertices comprising the face, shape = (*, 3)

    Returns
    -------
    np.ndarray
        center of the face, shape = (3,)

    """
    #mean unavailable in numba
    nverts = verts.shape[0]
    center = np.array((0., 0., 0.), dtype = 'float')
    center[0] = np.sum(verts[:,0])/nverts
    center[1] = np.sum(verts[:,1])/nverts
    center[2] = np.sum(verts[:,2])/nverts
    return center

#to avoid calling special scipy (binom) functions within numba
@numba.jit(nopython = True, cache = True, nogil = True)
def numberOfMonomials(dim:int,
                      order:int) -> int:
    """computes the number of monomials comprising the complete polynomial space of order order

    can be done with a simple scipy.special.binom(order + dim, dim)
    but unavailable from numba

    Parameters
    ----------
    dim : int
        dimension, number of variable comprising the polynomial space
    order : int
        poylnomial order

    Returns
    -------
    int
        number of monomials
    """
    if dim == 2:
        number = int( (order+2)*(order+1)/2 )
    elif dim ==3:
        number = int( (order+3)*(order+2)*(order+1)/6 )
    return number

@numba.jit(nopython=True, cache=True, nogil=True)
def evaluateMonomials(dim: int,
                      order: int,
                      verts: np.ndarray) -> np.ndarray:
    """evaluates all monomials comprising the complete polynomial space at vertices

    Parameters
    ----------
    dim : int
        dimension
    order : int
        polynomial order
    verts : np.ndarray
        coordiantes of vertices, shape (*, 3)

    Returns
    -------
    np.ndarray
        monomials evaluated at verts, shape (number of verts, number
        of monomials)

    """
    out = np.zeros((verts.shape[0],numberOfMonomials(dim, order)), dtype = 'float')
    restrict_complete = True
    ii = 0
    if dim == 2:
        for i in range(order+1):
            j_max = order-restrict_complete*i
            for j in range(j_max+1):
                out[:,ii] = verts[:,0]**i * verts[:,1]**j
                ii+=1
    else:
        for i in range(order+1):
            j_max = order-restrict_complete*i
            for j in range(j_max+1):
                k_max = order-restrict_complete*(i+j)
                for k in range(k_max+1):
                    out[:,ii] = verts[:,0]**i * verts[:,1]**j * verts[:,2]**k
                    ii+=1
    return out


@numba.jit(nopython = True, cache = True, nogil = True)
def getEdgeGeo(verts: np.ndarray,
               faceNormal: np.ndarray,
               faceCenter: np.ndarray) -> (np.ndarray, float, np.ndarray):
    """computes some geometric quantities associated with the edge

    computes the coordinates of the edge center, the distance from the
    edge center to the face center, the distance from the edge
    vertices to the edge center

    Parameters
    ----------
    verts : np.ndarray
        coordinates of vertices
    faceNormal : np.ndarray
        face normal, shape = (3,)
    faceCenter : np.ndarray
        face center, shape = (3,)

    Returns
    -------
    (np.ndarray, float, np.ndarray)

    edgeCenter, shape = (3,)
    edgeDistance, float
    vertsDistance, shape = (2,)
    """
    edgeNormal = np.cross(verts[1,:] - verts[0,:], faceNormal)
    edgeNormal = edgeNormal/np.sqrt(np.sum(edgeNormal**2))

    edgeCenter = (verts[0,:] + verts[1,:]) / 2
    edgeDistance = np.dot(edgeCenter - faceCenter, edgeNormal)

    vertsDistance = np.zeros(2, dtype = 'float')
    vertNormal = verts[0,:] - edgeCenter
    vertNormal = vertNormal/np.sqrt(np.sum(vertNormal**2))
    vertsDistance[0] =  vertNormal @ (verts[0,:] - edgeCenter).T
    vertNormal = verts[1,:] - edgeCenter
    vertNormal = vertNormal/np.sqrt(np.sum(vertNormal**2))
    vertsDistance[1] =  vertNormal @ (verts[1,:] - edgeCenter).T
    return edgeCenter, edgeDistance, vertsDistance

@numba.jit(nopython = True, cache = True, nogil = True)
def getMonomialIndex(dim:int,
                     order:int,
                     i:int,
                     j:int,
                     k:int) -> int:
    """get the index of a given monomial in the integral vector

    the integral is stored in a vector of shape (number of
    monomials,), this function helps you finding the index monomial
    x^i*y^j*z^k

    Parameters
    ----------
    dim : int
        dimension
    order : int
        polynomial order
    i : int
        power of x
    j : int
        power of y
    k : int
        power of z, k=0 if 2d

    Returns
    -------
    int
       index of the monomial    
    """
    if dim == 2:
        return getMonomialIndex(3, order, 0, i, j)
    ii = 0
    while (i>0):
        ii += (order+2)*(order+1)/2
        order-=1
        i-=1
    while (j>0):
        ii += order+1
        order-=1
        j-=1
    ii+=k
    return int(ii)

@numba.jit(nopython = True, cache = True, nogil = True)
def applyStokes(integral: np.ndarray,
                dim:int,
                order: int,
                edgeCenter: np.ndarray,
                geomDim: int):
    if dim == 2:
        ii = 0
        for i in range(order+1):
            for j in range(order+1-i):
                if i > 0:
                    integral[ii] += edgeCenter[0] * i * integral[getMonomialIndex(dim, order, i-1,  j, 0)]
                if j > 0:
                    integral[ii] += edgeCenter[1] * j * integral[getMonomialIndex(dim, order,  i, j-1, 0)]
                integral[ii] = integral[ii] / (geomDim+i+j)
                ii+=1
    if dim == 3:
        ii = 0
        for i in range(order+1):
            for j in range(order+1-i):
                for k in range(order + 1 - i - j):
                    if i > 0:
                        integral[ii] += edgeCenter[0] * i * integral[getMonomialIndex(dim, order, i-1,  j,   k)]
                    if j > 0:
                        integral[ii] += edgeCenter[1] * j * integral[getMonomialIndex(dim, order,  i, j-1,   k)]
                    if k > 0:
                        integral[ii] += edgeCenter[2] * k * integral[getMonomialIndex(dim, order,  i,   j, k-1)]
                    integral[ii] = integral[ii] / (geomDim+i+j+k)
                    ii+=1


@numba.jit(nopython = True, cache = True, nogil = True)
def integrateMonomialsEdge(dim: float,
                           order:float,
                           vertsPair:np.ndarray,
                           edgeCenter: np.ndarray,
                           edgeDistance: float,
                           vertsDistance: np.ndarray,
                           vertexPairMonomials: np.ndarray) -> np.ndarray:
    """integrate all monomial of the complete dim-vairate polynomial basis of order order on an edge 

    Parameters
    ----------
    dim : float
        dimension
    order : float
        polynomial order
    vertsPair : np.ndarray
        coordinates of vertices of the edge, shape (2,3)
    edgeCenter : np.ndarray
        coordinates of the center of the edge, shape (3,)
    edgeDistance : float
        distance from the edge to the center of the face
    vertsDistance : np.ndarray
        distance from vertices to the center of the edge, shape (2,)
    vertexPairMonomials : np.ndarray
        monomials evaluated at the vertices pair, shape (2, number of
        monomials)

    Returns
    -------
    np.ndarray
        integrated monomials vector, shape (number of monomials,)

    """
    integral  = vertsDistance[0]*vertexPairMonomials[0,:] + vertsDistance[1]*vertexPairMonomials[1,:]
    applyStokes(integral, dim, order, edgeCenter, 1)
    return integral

@numba.jit(nopython = True, cache = True, nogil = True)
def integrateMonomialsFace(dim: int,
                           order:int,
                           face: np.ndarray,
                           verts: np.ndarray,
                           vertexMonomials: np.ndarray) -> (np.ndarray, float):
    """integrate all monomials of the dim-variate complete polynomial basis or order order on a face

    Parameters
    ----------
    dim : int
        dimension
    order : int
        poylnomial order
    face : np.ndarray
        vector containing the index of vetices comprising the face,
        shape (number of edges of the face,)
    verts : np.ndarray
        coordintes of vertices of the whole poyltopal domain,
        shape (number of vertices, 3)
    vertexMonomials : np.ndarray
        monomials evaluated at all vertices of the polytopal domain,
        shape (number of verts, number of monomials)

    Returns
    -------
    (np.ndarray, float)
        integral, shape (number of monomials,)

    """
    numMono = vertexMonomials.shape[1]
    numEdges = face.size
    faceNormal = getFaceNormal(face, verts)
    faceCenter = getFaceCenter(verts[face,:])
    if dim == 2:
        faceDist = 0.
    if dim == 3:
        faceDist = np.dot(faceCenter, faceNormal)
    faceIntegral = np.zeros(numMono, dtype = 'float')
    for iEdge in range(numEdges):
        vertsPair = np.zeros((2,3), dtype = 'float')
        vertsPair[0,:] = verts[face[iEdge],:]
        vertsPair[1,:] = verts[face[(iEdge + 1) % numEdges],:]
        vertexPairMono = np.zeros((2,numMono), dtype = 'float')
        vertexPairMono[0,:] = vertexMonomials[face[iEdge],:]
        vertexPairMono[1,:] = vertexMonomials[face[(iEdge + 1) % numEdges],:]
        edgeCenter, edgeDistance, vertsDistance = getEdgeGeo(vertsPair, faceNormal, faceCenter)
        faceIntegral += edgeDistance * integrateMonomialsEdge(dim, order, vertsPair, edgeCenter, edgeDistance, vertsDistance, vertexPairMono)
    applyStokes(faceIntegral, dim, order, faceCenter, 2)
    return faceIntegral, faceDist

def integrateMonomialsPolygon(order: int,
                              face: np.ndarray,
                              verts: np.ndarray) -> np.ndarray:
    """computes the analytical integral of all monomials of the bi-variate poylnomial space of order order on a poylgon using Antonietti's algorithm

    IMPORTANT: assumes that points have already been mapped to the local bounding box. You can use /map_to_local_bb_2d/ for that.

    So far, there is no support for faces with holes, I yet have to come up with an idea to do it efficiently in python, sorry
    
    Parameters
    ----------
    order : int
        poylnomial order
    face : np.ndarray
        trigonomically orderd list of indeices of vertices comprising
        the polygon, shape (number of vertices,)
    verts : np.ndarray
        coordinates of vertices, shape (number of verts, 3)

    Returns
    -------
    np.ndarray
        integrated monomials, shape (number of monomials,)

    """
    vertexMonomials = evaluateMonomials(2, order, verts)
    integral, _ = integrateMonomialsFace(2, order, face, verts, vertexMonomials)
    return integral

@numba.jit(nopython = True, cache = True, nogil = True)
def delegate(order,dim,faces,verts,vertexMonomials):
    integral = np.zeros(int(numberOfMonomials(3, order)), dtype = 'float')
    for face in faces:
        faceIntegral,faceDist = integrateMonomialsFace(3, order,face, verts, vertexMonomials)
        integral += faceDist * faceIntegral

    ii = 0
    for i in range(order + 1):
        for j in range(order + 1 - i):
            for k in range(order + 1 - i - j):
                integral[ii] = integral[ii]/(3+i+j+k)
                ii+=1
    return integral

def integrateMonomialsPolyhedron(order: int,
                                 faces: list,
                                 verts: np.ndarray) -> np.ndarray:
    """computes the anlytical integral of all monomials of the tri-variate poylnomial space of order order on a polyhedron  using Antonietti's algorithm

    IMPORTANT: assumes that points have already been mapped to the local bounding box. You can use /map_to_local_bb_3d/ for that.

    Parameters
    ----------
    order : int
        polynomial order
    faces : list
        list of face, face contains trigonomically order indices of
        vertices comprising the face. Can be a list or an np.ndarray,
        np.ndarray is faster but only works if all faces have the same number of
        edges
    verts : np.ndarray
        coordinates of vertices, shape (number of vertices, 3)

    Returns
    -------
    np.ndarray
        integrated monomials, shape (number of monomials,)
    """
    vertexMonomials = evaluateMonomials(3, order, verts)
    integral = delegate(order, 3, faces, verts, vertexMonomials)
    return integral

if __name__=='__main__':
    # test using Paola's article as a reference
    #power of monomials in table 3
    monomials = [[5,5],[10,10],[20,20],[40,40],[10,5],[20,5],[40,5],[5,20],[5,40]]
    
    #p1 from Paola's article
    print('Testing on P1')
    face = np.array((0,1,2))
    verts = np.array([[-1.,1,-1],[-1,0,1],[0,0,0]]).T
    ref = [0,0.0111339078,0.0030396808,7.9534562047e-4,0,0,0,-0.005890191,-0.001868889]
    I = integrateMonomialsPolygon(80, face, verts)
    for ii,mono in enumerate(monomials):
        index_mono = getMonomialIndex(2,80,mono[0],mono[1],0)
        if abs(ref[ii]-I[index_mono])>1e-9:
            print(f'error too large on monomial {mono}')
