import cython
import numpy as np
cimport numpy as np

cdef extern from "P2Pwrap.h":

	ctypedef double REAL

	void direct_c(REAL *xi, int xiSize, REAL *yi, int yiSize, REAL *zi, int ziSize, REAL *xj, int xjSize, REAL *yj, int yjSize, REAL *zj, int zjSize, REAL *ds, int dsSize, REAL *dsx, int dsxSize, REAL *dsy, int dsySize, REAL *dsz, int dszSize, REAL *m, int mSize, REAL kappa, REAL eps)

	void P2P_c(REAL *xi, int xiSize, REAL *yi, int yiSize, REAL *zi, int ziSize, REAL *xj, int xjSize, REAL *yj, int yjSize, REAL *zj, int zjSize, int *interList, int interListSize, int *offTar, int offTarSize, int *sizeTar, int sizeTarSize, int *offSrc, int offSrcSize, int *offTwg, int offTwgSize, REAL *ds, int dsSize, REAL *dsx, int dsxSize, REAL *dsy, int dsySize, REAL *dsz, int dszSize, REAL *m, int mSize, REAL kappa, REAL eps)


def direct_c_cy(np.ndarray[REAL, ndim = 1, mode = "c"] xi, np.ndarray[REAL, ndim = 1, mode = "c"] yi, np.ndarray[REAL, ndim = 1, mode = "c"] zi, np.ndarray[REAL, ndim = 1, mode = "c"] xj, np.ndarray[REAL, ndim = 1, mode = "c"] yj, np.ndarray[REAL, ndim = 1, mode = "c"] zj, np.ndarray[REAL, ndim = 1, mode = "c"] ds, np.ndarray[REAL, ndim = 1, mode = "c"] dsx, np.ndarray[REAL, ndim = 1, mode = "c"] dsy, np.ndarray[REAL, ndim = 1, mode = "c"] dsz, np.ndarray[REAL, ndim = 1, mode = "c"] m, REAL kappa, REAL eps):
	cdef np.int32_t xiSize = len(xi)
	cdef np.int32_t yiSize = len(yi)
	cdef np.int32_t ziSize = len(zi)
	cdef np.int32_t xjSize = len(xj)
	cdef np.int32_t yjSize = len(yj)
	cdef np.int32_t zjSize = len(zj)
	cdef np.int32_t dsSize = len(ds)
	cdef np.int32_t dsxSize = len(dsx)
	cdef np.int32_t dsySize = len(dsy)
	cdef np.int32_t dszSize = len(dsz)
	cdef np.int32_t mSize = len(m)
	direct_c(<REAL*> &xi[0], <int> xiSize, <REAL*> &yi[0], <int> yiSize, <REAL*> &zi[0], <int> ziSize, <REAL*> &xj[0], <int> xjSize, <REAL*> &yj[0], <int> yjSize, <REAL*> &zj[0], <int> zjSize, <REAL*> &ds[0], <int> dsSize, <REAL*> &dsx[0], <int> dsxSize, <REAL*> &dsy[0], <int> dsySize, <REAL*> &dsz[0], <int> dszSize, <REAL*> &m[0], <int> mSize, <REAL> kappa, <REAL> eps)

def P2P_c_cy(np.ndarray[REAL, ndim = 1, mode = "c"] xi, np.ndarray[REAL, ndim = 1, mode = "c"] yi, np.ndarray[REAL, ndim = 1, mode = "c"] zi, np.ndarray[REAL, ndim = 1, mode = "c"] xj, np.ndarray[REAL, ndim = 1, mode = "c"] yj, np.ndarray[REAL, ndim = 1, mode = "c"] zj, np.ndarray[int, ndim = 1, mode = "c"] interList, np.ndarray[int, ndim = 1, mode = "c"] offTar, np.ndarray[int, ndim = 1, mode = "c"] sizeTar, np.ndarray[int, ndim = 1, mode = "c"] offSrc, np.ndarray[int, ndim = 1, mode = "c"] offTwg, np.ndarray[REAL, ndim = 1, mode = "c"] ds, np.ndarray[REAL, ndim = 1, mode = "c"] dsx, np.ndarray[REAL, ndim = 1, mode = "c"] dsy, np.ndarray[REAL, ndim = 1, mode = "c"] dsz, np.ndarray[REAL, ndim = 1, mode = "c"] m, REAL kappa, REAL eps):
	cdef np.int32_t xiSize = len(xi)
	cdef np.int32_t yiSize = len(yi)
	cdef np.int32_t ziSize = len(zi)
	cdef np.int32_t xjSize = len(xj)
	cdef np.int32_t yjSize = len(yj)
	cdef np.int32_t zjSize = len(zj)
	cdef np.int32_t dsSize = len(ds)
	cdef np.int32_t dsxSize = len(dsx)
	cdef np.int32_t dsySize = len(dsy)
	cdef np.int32_t dszSize = len(dsz)
	cdef np.int32_t mSize = len(m)
	cdef np.int32_t interListSize = len(interList)
	cdef np.int32_t offTarSize = len(offTar)
	cdef np.int32_t sizeTarSize = len(sizeTar)
	cdef np.int32_t offSrcSize = len(offSrc)
	cdef np.int32_t offTwgSize = len(offTwg)
	P2P_c(<REAL*> &xi[0], <int> xiSize, <REAL*> &yi[0], <int> yiSize, <REAL*> &zi[0], <int> ziSize, <REAL*> &xj[0], <int> xjSize, <REAL*> &yj[0], <int> yjSize, <REAL*> &zj[0], <int> zjSize, <int*> &interList[0], <int> interListSize, <int*> &offTar[0], <int> offTarSize, <int*> &sizeTar[0], <int> sizeTarSize, <int*> &offSrc[0], <int> offSrcSize, <int*> &offTwg[0], <int> offTwgSize, <REAL*> &ds[0], <int> dsSize, <REAL*> &dsx[0], <int> dsxSize, <REAL*> &dsy[0], <int> dsySize, <REAL*> &dsz[0], <int> dszSize, <REAL*> &m[0], <int> mSize, <REAL> kappa, <REAL> eps)
