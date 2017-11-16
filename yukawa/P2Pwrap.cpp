#include <stdio.h>
#include <time.h>
#include <cmath>
#include <cstdlib>
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
              REAL kappa, REAL eps)
{   
    /*
    double kappa = 1.5;

    double eps = 1e-10;
    double *xi = new double[N];
    double *yi = new double[N];
    double *zi = new double[N];
    double *xj = new double[N];
    double *yj = new double[N];
    double *zj = new double[N];
    double *m = new double[N];
    double *ds = new double[N];
    double *dsx = new double[N];
    double *dsy = new double[N];
    double *dsz = new double[N];

    for (int i; i<N; i++)
    {
        xi[i] = rand()/(1.+RAND_MAX);
        yi[i] = rand()/(1.+RAND_MAX);
        zi[i] = rand()/(1.+RAND_MAX);
        xj[i] = rand()/(1.+RAND_MAX);
        yj[i] = rand()/(1.+RAND_MAX);
        zj[i] = rand()/(1.+RAND_MAX);
        m[i] = 1.0f/N;
    }   
    clock_t begin = clock();
    */

    // Direct summation
    REAL r, r2, dx, dy, dz; 
    //time_t begin, end;
    //time (&begin);
    for (int i=0; i<xiSize; i++)
    {   
        //ds[i] = -m[i]/sqrt(eps);
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
    /*
    clock_t end = clock();
    //time (&end);
    //double time_sec = difftime(end,begin);

    double time_sec = double(end-begin)/CLOCKS_PER_SEC;
//    time_sec /= 100;

    printf ("Time: %fs\n",time_sec);
    return 0;
    */
}


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
