%module M2Pwrap

%{
#define SWIG_FILE_WITH_INIT
extern void M2P_c(double *p, int pSize, 
                    double *px, int pxSize, 
                    double *py, int pySize, 
                    double *pz, int pzSize, 
                    double *mult, int multSize, 
                    double *xa, int xaSize,
                    double *ya, int yaSize, 
                    double *za, int zaSize, 
                    double *xc_arr, int xcSize,
                    double *yc_arr, int ycSize, 
                    double *zc_arr, int zcSize, 
                    int *offMlt, int offMltSize,
                    int *offTar, int offTarSize,
                    int *sizeTar, int sizeTarSize,
                    int P, double kappa, int Nm, int *Index, int indSize);
extern void getCoeff(double *a, int aSize, 
                    double *ax, int axSize, 
                    double *ay, int aySize, 
                    double *az, int azSize, 
                    double dx, double dy, double dz, int P, double kappa, int *Index, int indSize);
extern void setIndex_arr(int P, int N, 
                        int *indices, int indicesSize,
                        int * ii    , int iiSize,
                        int * jj    , int jjSize,
                        int * kk    , int kkSize);
extern int setIndex(int P, int i, int j, int k);
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double* INPLACE_ARRAY1, int DIM1){(double *a, int aSize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *ax, int axSize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *ay, int aySize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *az, int azSize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *p, int pSize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *px, int pxSize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *py, int pySize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *pz, int pzSize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *mult, int multSize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *xa, int xaSize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *ya, int yaSize)};
%apply (double* INPLACE_ARRAY1, int DIM1){(double *za, int zaSize)};
%apply (double* IN_ARRAY1, int DIM1){(double *xc_arr, int xcSize)};
%apply (double* IN_ARRAY1, int DIM1){(double *yc_arr, int ycSize)};
%apply (double* IN_ARRAY1, int DIM1){(double *zc_arr, int zcSize)};
%apply (int* IN_ARRAY1, int DIM1){(int *indices, int indicesSize)};
%apply (int* IN_ARRAY1, int DIM1){(int *ii, int iiSize)};
%apply (int* IN_ARRAY1, int DIM1){(int *jj, int jjSize)};
%apply (int* IN_ARRAY1, int DIM1){(int *kk, int kkSize)};
%apply (int* IN_ARRAY1, int DIM1){(int *offMlt, int offMltSize)};
%apply (int* IN_ARRAY1, int DIM1){(int *offTar, int offTarSize)};
%apply (int* IN_ARRAY1, int DIM1){(int *sizeTar, int sizeTarSize)};
%apply (int* IN_ARRAY1, int DIM1){(int *Index, int indSize)};
%apply (int* IN_ARRAY1, int DIM1){(int *target, int tarSize)};
extern void M2P_c(double *p, int pSize, 
                    double *px, int pxSize, 
                    double *py, int pySize, 
                    double *pz, int pzSize, 
                    double *mult, int multSize, 
                    double *xa, int xaSize,
                    double *ya, int yaSize, 
                    double *za, int zaSize, 
                    double *xc_arr, int xcSize,
                    double *yc_arr, int ycSize, 
                    double *zc_arr, int zcSize, 
                    int *offMlt, int offMltSize,
                    int *offTar, int offTarSize,
                    int *sizeTar, int sizeTarSize,
                    int P, double kappa, int Nm, int *Index, int indSize);
extern void getCoeff(double *a, int aSize, 
                    double *ax, int axSize, 
                    double *ay, int aySize, 
                    double *az, int azSize, 
                    double dx, double dy, double dz, int P, double kappa, int *Index, int indSize);
extern void setIndex_arr(int P, int N, 
                        int *indices, int indicesSize,
                        int * ii    , int iiSize,
                        int * jj    , int jjSize,
                        int * kk    , int kkSize);
extern int setIndex(int P, int i, int j, int k);
%clear (double *a, int aSize);
%clear (double *ax, int axSize);
%clear (double *ay, int aySize);
%clear (double *az, int azSize);
%clear (double *p, int pSize);
%clear (double *px, int pxSize);
%clear (double *py, int pySize);
%clear (double *pz, int pzSize);
%clear (double *mult, int multSize);
%clear (double *xa, int xaSize);
%clear (double *ya, int yaSize);
%clear (double *za, int zaSize);
%clear (double *xc_arr, int xcSize);
%clear (double *yc_arr, int ycSize);
%clear (double *zc_arr, int zcSize);
%clear (int *indices, int indicesSize);
%clear (int *ii, int iiSize);
%clear (int *jj, int jjSize);
%clear (int *kk, int kkSize);
%clear (int *offMlt, int offMltSize);
%clear (int *offTar, int offTarSize);
%clear (int *sizeTar, int sizeTarSize);
%clear (int *Index, int indSize);
%clear (int *target, int tarSize);
