#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#define REAL double

void direct_c(REAL *xi, int xiSize, 
              REAL *yi, int yiSize, 
              REAL *zi, int ziSize, 
              REAL *xj, int xjSize, 
              REAL *yj, int yjSize, 
              REAL *zj, int zjSize, 
              REAL *ds, int dsSize, 
              REAL *dsx, int dsxSize, 
              REAL *dsy, int dsySize, 
              REAL *dsz, int dszSize, 
              REAL *m, int mSize, 
              REAL kappa, REAL eps);

void P2P_c(REAL *xi, int xiSize, 
              REAL *yi, int yiSize, 
              REAL *zi, int ziSize, 
              REAL *xj, int xjSize, 
              REAL *yj, int yjSize, 
              REAL *zj, int zjSize, 
              int *interList, int interListSize,
              int *offTar, int offTarSize,
              int *sizeTar, int sizeTarSize,
              int *offSrc, int offSrcSize,
              int *offTwg, int offTwgSize, 
              REAL *ds, int dsSize, 
              REAL *dsx, int dsxSize, 
              REAL *dsy, int dsySize, 
              REAL *dsz, int dszSize, 
              REAL *m, int mSize, 
              REAL kappa, REAL eps);

void direct_c(REAL *xi, int xiSize, 
              REAL *yi, int yiSize, 
              REAL *zi, int ziSize, 
              REAL *xj, int xjSize, 
              REAL *yj, int yjSize, 
              REAL *zj, int zjSize, 
              REAL *ds, int dsSize, 
              REAL *dsx, int dsxSize, 
              REAL *dsy, int dsySize, 
              REAL *dsz, int dszSize, 
              REAL *m, int mSize, 
              REAL kappa, REAL eps)
{   
    REAL r, r2, dx, dy, dz; 

    for (int i=0; i<xiSize; i++)
    {   
        for (int j=0; j<xjSize; j++)
        {   
            dx = xi[i] - xj[j]; 
            dy = yi[i] - yj[j]; 
            dz = zi[i] - zj[j]; 
            r = sqrt(dx*dx + dy*dy + dz*dz + eps);
            r2 = r*r;
            ds[i] += m[j]*exp(-kappa*r)/r;
            dsx[i] += -m[j]*dx*exp(-kappa*r)/r2*(kappa+1/r);
            dsy[i] += -m[j]*dy*exp(-kappa*r)/r2*(kappa+1/r);
            dsz[i] += -m[j]*dz*exp(-kappa*r)/r2*(kappa+1/r);

        }   
    }   
};


void P2P_c(REAL *xi, int xiSize, 
              REAL *yi, int yiSize, 
              REAL *zi, int ziSize, 
              REAL *xj, int xjSize, 
              REAL *yj, int yjSize, 
              REAL *zj, int zjSize, 
              int *interList, int interListSize,
              int *offTar, int offTarSize,
              int *sizeTar, int sizeTarSize,
              int *offSrc, int offSrcSize,
              int *offTwg, int offTwgSize, 
              REAL *ds, int dsSize, 
              REAL *dsx, int dsxSize, 
              REAL *dsy, int dsySize, 
              REAL *dsz, int dszSize, 
              REAL *m, int mSize, 
              REAL kappa, REAL eps)
{  
    REAL r, r2, dx, dy, dz, aux, sum_p, sum_x, sum_y, sum_z, x, y, z; 
    int CI_start, CI_end, CJ_start, CJ_end, list_start, list_end, CJ;

    for (int tarTwg=0; tarTwg<offTarSize; tarTwg++)
    {   
        CI_start = offTar[tarTwg];
        CI_end   = CI_start + sizeTar[tarTwg];
        list_start = offTwg[tarTwg];
        list_end   = offTwg[tarTwg+1];

        for (int i=CI_start; i<CI_end; i++)
        {
            sum_p = 0.;
            sum_x = 0.;
            sum_y = 0.;
            sum_z = 0.;
            x = xi[i];
            y = yi[i];
            z = zi[i];

            for (int lst=list_start; lst<list_end; lst++)
            {
                CJ = interList[lst];
                CJ_start = offSrc[CJ];
                CJ_end = offSrc[CJ+1];

                for (int j=CJ_start; j<CJ_end; j++)
                {   
                    dx = x - xj[j]; 
                    dy = y - yj[j]; 
                    dz = z - zj[j]; 
                    r2 = dx*dx + dy*dy + dz*dz + eps;
                    r  = sqrt(r2);
                    aux     = m[j]*exp(-kappa*r)/r;
                    sum_p  += aux;

                    aux *= (kappa*r+1)/r2;
                    sum_x += -dx*aux;
                    sum_y += -dy*aux;
                    sum_z += -dz*aux;
                }
            }

            ds[i] += sum_p;
            dsx[i] += sum_x;
            dsy[i] += sum_y;
            dsz[i] += sum_z;
        }
    }   
}
