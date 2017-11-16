import cython
import numpy as np
cimport numpy as np

cdef extern from "M2Pwrap.h":

	ctypedef double REAL

	int setIndex(int P, int i, int j, int k)

	void setIndex_arr(int P, int N, int *indices, int indicesSize, int *ii, int iiSize, int *jj, int jjSize, int *kk, int kkSize)

	void getCoeff(REAL *a, int aSize, REAL *ax, int axSize, REAL *ay, int aySize, REAL *az, int azSize, REAL dx, REAL dy, REAL dz, int P, REAL kappa, int *Index, int indSize)

	void M2P_c(REAL *p, int pSize, REAL *px, int pxSize, REAL *py, int pySize, REAL *pz, int pzSize, REAL *mult, int multSize, REAL *x, int xSize, REAL *y, int ySize, REAL *z, int zSize, REAL *xc, int xcSize, REAL *yc, int ycSize, REAL *zc, int zcSize, int *offMlt, int offMltSize, int *offTar, int offTarSize, int *sizeTar, int sizeTarSize, int P, REAL kappa, int Nm, int *Index, int indSize)


def setIndex_cy(int P, int i, int j, int k):
	return setIndex(<int> P, <int> i, <int> j, <int> k)

def setIndex_arr_cy(int P, int N, np.ndarray[int, ndim = 1, mode = "c"] indices, np.ndarray[int, ndim = 1, mode = "c"] ii, np.ndarray[int, ndim = 1, mode = "c"] jj, np.ndarray[int, ndim = 1, mode = "c"] kk):
	cdef np.int32_t indicesSize = len(indices)
	cdef np.int32_t iiSize = len(ii)
	cdef np.int32_t jjSize = len(jj)
	cdef np.int32_t kkSize = len(kk)
	setIndex_arr(<int> P, <int> N, <int*> &indices[0], <int> indicesSize, <int*> &ii[0], <int> iiSize, <int*> &jj[0], <int> jjSize, <int*> &kk[0], <int> kkSize)

def getCoeff_cy(np.ndarray[REAL, ndim = 1, mode = "c"] a, np.ndarray[REAL, ndim = 1, mode = "c"] ax, np.ndarray[REAL, ndim = 1, mode = "c"] ay, np.ndarray[REAL, ndim = 1, mode = "c"] az, REAL dx, REAL dy, REAL dz, int P, REAL kappa, np.ndarray[long, ndim = 1, mode = "c"] Index):
	cdef np.int32_t aSize = len(a)
	cdef np.int32_t axSize = len(ax)
	cdef np.int32_t aySize = len(ay)
	cdef np.int32_t azSize = len(az)
	cdef np.int32_t indSize = len(Index)
	getCoeff(<REAL*> &a[0], <int> aSize, <REAL*> &ax[0], <int> axSize, <REAL*> &ay[0], <int> aySize, <REAL*> &az[0], <int> azSize, <REAL> dx, <REAL> dy, <REAL> dz, <int> P, <REAL> kappa, <int*> &Index[0], <int> indSize)

def M2P_c_cy(np.ndarray[REAL, ndim = 1, mode = "c"] p, np.ndarray[REAL, ndim = 1, mode = "c"] px, np.ndarray[REAL, ndim = 1, mode = "c"] py, np.ndarray[REAL, ndim = 1, mode = "c"] pz, np.ndarray[REAL, ndim = 1, mode = "c"] mult, np.ndarray[REAL, ndim = 1, mode = "c"] x, np.ndarray[REAL, ndim = 1, mode = "c"] y, np.ndarray[REAL, ndim = 1, mode = "c"] z, np.ndarray[REAL, ndim = 1, mode = "c"] xc, np.ndarray[REAL, ndim = 1, mode = "c"] yc, np.ndarray[REAL, ndim = 1, mode = "c"] zc, np.ndarray[int, ndim = 1, mode = "c"] offMlt, np.ndarray[int, ndim = 1, mode = "c"] offTar, np.ndarray[int, ndim = 1, mode = "c"] sizeTar, int P, REAL kappa, int Nm, np.ndarray[int, ndim = 1, mode = "c"] Index):
	cdef np.int32_t pSize = len(p)
	cdef np.int32_t pxSize = len(px)
	cdef np.int32_t pySize = len(py)
	cdef np.int32_t pzSize = len(pz)
	cdef np.int32_t multSize = len(mult)
	cdef np.int32_t xSize = len(x)
	cdef np.int32_t ySize = len(y)
	cdef np.int32_t zSize = len(z)
	cdef np.int32_t xcSize = len(xc)
	cdef np.int32_t ycSize = len(yc)
	cdef np.int32_t zcSize = len(zc)
	cdef np.int32_t offMltSize = len(offMlt)
	cdef np.int32_t offTarSize = len(offTar)
	cdef np.int32_t sizeTarSize = len(sizeTar)
	cdef np.int32_t indSize = len(Index)
	M2P_c(<REAL*> &p[0], <int> pSize, <REAL*> &px[0], <int> pxSize, <REAL*> &py[0], <int> pySize, <REAL*> &pz[0], <int> pzSize, <REAL*> &mult[0], <int> multSize, <REAL*> &x[0], <int> xSize, <REAL*> &y[0], <int> ySize, <REAL*> &z[0], <int> zSize, <REAL*> &xc[0], <int> xcSize, <REAL*> &yc[0], <int> ycSize, <REAL*> &zc[0], <int> zcSize, <int*> &offMlt[0], <int> offMltSize, <int*> &offTar[0], <int> offTarSize, <int*> &sizeTar[0], <int> sizeTarSize, <int> P, <REAL> kappa, <int> Nm, <int*> &Index[0], <int> indSize)
