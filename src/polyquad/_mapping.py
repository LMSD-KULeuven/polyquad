import numpy as np

def map_to_local_bb_2d(verts: np.ndarray) -> (np.ndarray, float):
    """map the element to the reference bounding box and computes the jacobian

    Parameters
    ----------
    verts : np.ndarray
        coordinates of vertices, shape = ( *, 2)

    Returns
    -------
    (np.ndarray, float)
        (mapped coordinates, jacobian)

    """
    #span of the element
    dx = abs(verts[:,0].max() - verts[:,0].min())
    dy = abs(verts[:,1].max() - verts[:,1].min())
    #stretch matrix
    S = np.array(((2/dx,    0),
                  (   0, 2/dy)))
    #perform strech
    v = (S @ verts.T).T
    #perform transltation
    translate = np.min(v, axis = 0) + 1
    v -= translate
    #jacobian
    j = dx * dy
    return v, j

def map_to_local_bb_3d(verts: np.ndarray) -> (np.ndarray, float):
    """map the element to the reference bounding box and computes the jacobian

    Parameters
    ----------
    verts : np.ndarray
        coordinates of vertices, shape = ( *, 3)

    Returns
    -------
    (np.ndarray, float)
        (mapped coordinates, jacobian)

    """
    #span of the element
    dx = abs(verts[:,0].max() - verts[:,0].min())
    dy = abs(verts[:,1].max() - verts[:,1].min())
    dz = abs(verts[:,2].max() - verts[:,2].min())
    #stretch matrix
    S = np.array(((2/dx,    0,    0),
                  (   0, 2/dy,    0),
                  (   0,    0, 2/dz)))
    #perform strech
    v = (S @ verts.T).T
    #perform transltation
    translate = np.min(v, axis = 0) + 1
    v -= translate
    #jacobian
    j = dx * dy * dz
    return v, j
