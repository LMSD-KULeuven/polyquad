import numpy as np
import scipy as sp

QRPV={}

def mkDecomposition(order: int):
    """generates points, conputes the vandermonde matrix, computes QR decomposition and store everything

    Parameters
    ----------
    order : int
        polynomial order of the integrand

    """
    jj,_ = sp.special.roots_legendre(order+1)
    X,Y,Z = np.meshgrid(jj,jj,jj)
    x,y,z = X.flatten(), Y.flatten(), Z.flatten()
    x,y,z = x[np.newaxis],y[np.newaxis],z[np.newaxis]
    xpow,ypow,zpow = trios(order)
    vt = x**xpow * y**ypow * z**zpow
    q,r,p = sp.linalg.qr(vt, mode = 'economic', pivoting = True)
    M,N = r.shape
    RinvQt = np.linalg.inv(r[:M,:M]) @ q.T
    pts = np.hstack((x.T,y.T,z.T))
    pts = pts[:M,:]
    QRPV[order] = [RinvQt, p, vt.T, pts]

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
    xPow, yPow, zPow = v1[mask], v2[mask], v3[mask]
    return yPow.reshape((yPow.size,1)), xPow.reshape((xPow.size,1)),zPow.reshape((zPow.size,1))

# def reshape_monomials(mono,order):
#     xp, yp, zp = trios(order)
#     m = np.zeros((len(xp),1),dtype = 'float')
#     for ii in range(len(xp)):
#         m[ii,0] = mono[xp[ii],yp[ii],zp[ii]]
#     return m

def moment_matching(order: int,
                    mono: np.ndarray,
                    residual = False) -> (np.ndarray,np.ndarray):
    """perform moment matching with analytical integral of monomials

    Parameters
    ----------
    order : int
        polynomial order of the integrand
    mono : np.ndarray
        integrated monomials, shape = (TO BE DEFINED)
    residual : bool
        returns the the residual if true

    Returns
    -------
    pts : np.ndarray
      coordinates of integration points
    w : np.ndarray
      weights associated with integraiton points
    residual: float
      quality indicator of the quadrature, see equation 16
    """
    if not(order in QRPV.keys()):
        print('decomposing')
        mkDecomposition(order)
    #extract decomposition
    RinvQt,p,v,pts = QRPV[order]
    N,M = v.shape
    m = mono
    y = RinvQt @ mono
    if residual:
        x = np.zeros(N, dtype = 'float')
        x[p[:M]] = y
        res = np.linalg.norm(v.T @ x - m,2)
        return pts, y, res
    return pts,y
