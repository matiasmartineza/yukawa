import cython
import numpy as np
cimport numpy as np

cdef extern from "calculateMultipoles.h":

	ctypedef double REAL

	void P2M(REAL *M, int Msize, REAL *x, int xSize, REAL *y, int ySize, REAL *z, int zSize, REAL *m, int mSize, REAL xc, REAL yc, REAL zc, int *I, int Isize, int *J, int Jsize, int *K, int Ksize)

	void M2M(REAL *MP, int MPsize, REAL *MC, int MCsize, REAL dx, REAL dy, REAL dz, int *I, int Isize, int *J, int Jsize, int *K, int Ksize, REAL *cI, int cIsize, REAL *cJ, int cJsize, REAL *cK, int cKsize, int *Imi, int Imisize, int *Jmj, int Jmjsize, int *Kmk, int Kmksize, int *index, int indexSize, int *ptr, int ptrSize)

def P2M_cy(np.ndarray[REAL, ndim = 1, mode = "c"] M, np.ndarray[REAL, ndim = 1, mode = "c"] x, np.ndarray[REAL, ndim = 1, mode = "c"] y, np.ndarray[REAL, ndim = 1, mode = "c"] z, np.ndarray[REAL, ndim = 1, mode = "c"] m, REAL xc, REAL yc, REAL zc, np.ndarray[int, ndim = 1, mode = "c"] I, np.ndarray[int, ndim = 1, mode = "c"] J, np.ndarray[int, ndim = 1, mode = "c"] K):
	cdef np.int32_t Msize = len(M)
	cdef np.int32_t xSize = len(x)
	cdef np.int32_t ySize = len(y)
	cdef np.int32_t zSize = len(z)
	cdef np.int32_t mSize = len(m)
	cdef np.int32_t Isize = len(I)
	cdef np.int32_t Jsize = len(J)
	cdef np.int32_t Ksize = len(K)
	P2M(<REAL*> &M[0], <int> Msize, <REAL*> &x[0], <int> xSize, <REAL*> &y[0], <int> ySize, <REAL*> &z[0], <int> zSize, <REAL*> &m[0], <int> mSize, <REAL> xc, <REAL> yc, <REAL> zc, <int*> &I[0], <int> Isize, <int*> &J[0], <int> Jsize, <int*> &K[0], <int> Ksize)

def M2M_cy(np.ndarray[REAL, ndim = 1, mode = "c"] MP, np.ndarray[REAL, ndim = 1, mode = "c"] MC, REAL dx, REAL dy, REAL dz, np.ndarray[int, ndim = 1, mode = "c"] I, np.ndarray[int, ndim = 1, mode = "c"] J, np.ndarray[int, ndim = 1, mode = "c"] K, np.ndarray[REAL, ndim = 1, mode = "c"] cI, np.ndarray[REAL, ndim = 1, mode = "c"] cJ, np.ndarray[REAL, ndim = 1, mode = "c"] cK, np.ndarray[int, ndim = 1, mode = "c"] Imi, np.ndarray[int, ndim = 1, mode = "c"] Jmj, np.ndarray[int, ndim = 1, mode = "c"] Kmk, np.ndarray[int, ndim = 1, mode = "c"] index, np.ndarray[int, ndim = 1, mode = "c"] ptr):
	cdef np.int32_t MPsize = len(MP)
	cdef np.int32_t MCsize = len(MC)
	cdef np.int32_t Isize = len(I)
	cdef np.int32_t Jsize = len(J)
	cdef np.int32_t Ksize = len(K)
	cdef np.int32_t cIsize = len(cI)
	cdef np.int32_t cJsize = len(cJ)
	cdef np.int32_t cKsize = len(cK)
	cdef np.int32_t Imisize = len(Imi)
	cdef np.int32_t Jmjsize = len(Jmj)
	cdef np.int32_t Kmksize = len(Kmk)
	cdef np.int32_t indexSize = len(index)
	cdef np.int32_t ptrSize = len(ptr)
	M2M(<REAL*> &MP[0], <int> MPsize, <REAL*> &MC[0], <int> MCsize, <REAL> dx, <REAL> dy, <REAL> dz, <int*> &I[0], <int> Isize, <int*> &J[0], <int> Jsize, <int*> &K[0], <int> Ksize, <REAL*> &cI[0], <int> cIsize, <REAL*> &cJ[0], <int> cJsize, <REAL*> &cK[0], <int> cKsize, <int*> &Imi[0], <int> Imisize, <int*> &Jmj[0], <int> Jmjsize, <int*> &Kmk[0], <int> Kmksize, <int*> &index[0], <int> indexSize, <int*> &ptr[0], <int> ptrSize)
