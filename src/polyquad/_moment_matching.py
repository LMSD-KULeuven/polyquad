import numpy as np
import scipy as sp

class VandermondeCollection:
    def __init__(self):
        self.order_dict = {}
    def done(self, order):
        return order in self.order_dict.keys()
    def add(self, order):
        self.order_dict[order] = Vandermonde(order)

class Vandermonde:
    def __init__(self, order):
        self.mkVandermonde(order)
        self.decompose()
    def mkVandermonde(self, order):
        jj,_ = sp.special.roots_legendre(order+1)
        X,Y,Z = np.meshgrid(jj,jj,jj)
        x,y,z = X.flatten(), Y.flatten(), Z.flatten()
        x,y,z = x[np.newaxis],y[np.newaxis],z[np.newaxis]
        xpow,ypow,zpow = trios(order)
        vt = x**xpow * y**ypow * z**zpow
        self.matrice = vt
        self.pts = np.hstack((x.T,y.T,z.T))
    def decompose(self):
        q,r,p = sp.linalg.qr(self.matrice, mode = 'economic', pivoting = True)
        M, _ = self.matrice.shape
        self.q = q
        # self.r = r
        self.p = p
        # self.r_truncated = r[:M,:M]
        self.pts_truncated = self.pts[p[:M],:]
        lu, piv = sp.linalg.lu_factor(r[:M,:M])
        self.lu, self.piv = lu,piv

collec = VandermondeCollection()

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
    if not(collec.done(order)):
        collec.add(order)

    lu,piv = collec.order_dict[order].lu, collec.order_dict[order].piv
    q = collec.order_dict[order].q
    y = sp.linalg.lu_solve((lu,piv), q.T @ mono)
    if residual:
        M,N = collec.order_dict[order].matrice.shape
        x = np.zeros(N, dtype = 'float')
        p = collec.order_dict[order].p
        x[p[:M]] = y
        res = np.linalg.norm(collec.order_dict[order].matrice @ x - mono,2)/ np.linalg.norm(mono, 2)
        return collec.order_dict[order].pts_truncated, y, res
    return collec.order_dict[order].pts_truncated,y
