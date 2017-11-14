import cython
import numpy as np
cimport numpy as np

cdef extern from "P2Pwrap.h":

	ctypedef double REAL

	void direct_c(REAL *xi, int xiSize, REAL *yi, int yiSize, REAL *zi, int ziSize, REAL *xj, int xjSize, REAL *yj, int yjSize, REAL *zj, int zjSize, REAL *ds, int dsSize, REAL *dsx, int dsxSize, REAL *dsy, int dsySize, REAL *dsz, int dszSize, REAL *m, int mSize, REAL kappa, REAL eps)

	void P2P_c(REAL *xi, int xiSize, REAL *yi, int yiSize, REAL *zi, int ziSize, REAL *xj, int xjSize, REAL *yj, int yjSize, REAL *zj, int zjSize, int *interList, int interListSize, int *offTar, int offTarSize, int *sizeTar, int sizeTarSize, int *offSrc, int offSrcSize, int *offTwg, int offTwgSize, REAL *ds, int dsSize, REAL *dsx, int dsxSize, REAL *dsy, int dsySize, REAL *dsz, int dszSize, REAL *m, int mSize, REAL kappa, REAL eps)


def direct_c_cy(np.ndarray[REAL, ndim = 1, mode = "c"] xi, np.ndarray[REAL, ndim = 1, mode = "c"] yi, np.ndarray[REAL, ndim = 1, mode = "c"] zi, np.ndarray[REAL, ndim = 1, mode = "c"] xj, np.ndarray[REAL, ndim = 1, mode = "c"] yj, np.ndarray[REAL, ndim = 1, mode = "c"] zj, np.ndarray[REAL, ndim = 1, mode = "c"] ds, np.ndarray[REAL, ndim = 1, mode = "c"] dsx, np.ndarray[REAL, ndim = 1, mode = "c"] dsy, np.ndarray[REAL, ndim = 1, mode = "c"] dsz, np.ndarray[REAL, ndim = 1, mode = "c"] m, REAL kappa, REAL eps):
	cdef np.int32_t n = len(xi)
	direct_c(<REAL*> &xi[0], <int> n, <REAL*> &yi[0], <int> n, <REAL*> &zi[0], <int> n, <REAL*> &xj[0], <int> n, <REAL*> &yj[0], <int> n, <REAL*> &zj[0], <int> n, <REAL*> &ds[0], <int> n, <REAL*> &dsx[0], <int> n, <REAL*> &dsy[0], <int> n, <REAL*> &dsz[0], <int> n, <REAL*> &m[0], <int> n, <REAL> kappa, <REAL> eps)

def P2P_c_cy(np.ndarray[REAL, ndim = 1, mode = "c"] xi, np.ndarray[REAL, ndim = 1, mode = "c"] yi, np.ndarray[REAL, ndim = 1, mode = "c"] zi, np.ndarray[REAL, ndim = 1, mode = "c"] xj, np.ndarray[REAL, ndim = 1, mode = "c"] yj, np.ndarray[REAL, ndim = 1, mode = "c"] zj, np.ndarray[int, ndim = 1, mode = "c"] interList, np.ndarray[int, ndim = 1, mode = "c"] offTar, np.ndarray[int, ndim = 1, mode = "c"] sizeTar, np.ndarray[int, ndim = 1, mode = "c"] offSrc, np.ndarray[int, ndim = 1, mode = "c"] offTwg, np.ndarray[REAL, ndim = 1, mode = "c"] ds, np.ndarray[REAL, ndim = 1, mode = "c"] dsx, np.ndarray[REAL, ndim = 1, mode = "c"] dsy, np.ndarray[REAL, ndim = 1, mode = "c"] dsz, np.ndarray[REAL, ndim = 1, mode = "c"] m, REAL kappa, REAL eps):
	cdef np.int32_t n = len(xi)
	cdef np.int32_t interListSize = len(interList)
	cdef np.int32_t offTarSize = len(offTar)
	cdef np.int32_t sizeTarSize = len(sizeTar)
	cdef np.int32_t offSrcSize = len(offSrc)
	cdef np.int32_t offTwgSize = len(offTwg)
	P2P_c(<REAL*> &xi[0], <int> n, <REAL*> &yi[0], <int> n, <REAL*> &zi[0], <int> n, <REAL*> &xj[0], <int> n, <REAL*> &yj[0], <int> n, <REAL*> &zj[0], <int> n, <int*> &interList[0], <int> interListSize, <int*> &offTar[0], <int> offTarSize, <int*> &sizeTar[0], <int> sizeTarSize, <int*> &offSrc[0], <int> offSrcSize, <int*> &offTwg[0], <int> offTwgSize, <REAL*> &ds[0], <int> n, <REAL*> &dsx[0], <int> n, <REAL*> &dsy[0], <int> n, <REAL*> &dsz[0], <int> n, <REAL*> &m[0], <int> n, <REAL> kappa, <REAL> eps)
