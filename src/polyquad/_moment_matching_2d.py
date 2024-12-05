import numpy as np
import scipy as sp

from ._utils import duos

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
        X,Y = np.meshgrid(jj,jj)
        x,y = X.flatten(), Y.flatten()
        x,y = x[np.newaxis],y[np.newaxis]#,z[np.newaxis]
        xpow,ypow = duos(order)
        vt = x**xpow * y**ypow
        self.matrice = vt
        self.pts = np.hstack((x.T,y.T))
    def decompose(self):
        q,r,p = sp.linalg.qr(self.matrice, mode = 'economic', pivoting = True)
        M, _ = self.matrice.shape
        self.q = q
        self.p = p
        self.r_truncated = np.asfortranarray(r[:M,:M]) # converting to fortran speeds up the backward substitution step
        self.pts_truncated = self.pts[p[:M],:]

collec = VandermondeCollection()


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

    r = collec.order_dict[order].r_truncated
    q = collec.order_dict[order].q
    y,_ = sp.linalg.lapack.dtrtrs(r, q.T @ mono)
    if residual:
        M,N = collec.order_dict[order].matrice.shape
        x = np.zeros(N, dtype = 'float')
        p = collec.order_dict[order].p
        x[p[:M]] = y
        res = np.linalg.norm(collec.order_dict[order].matrice @ x - mono,2)/ np.linalg.norm(mono, 2)
        return collec.order_dict[order].pts_truncated, y, res
    return collec.order_dict[order].pts_truncated,y
