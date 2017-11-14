#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#define REAL double

double get_time (void);

int setIndex(int P, int i, int j, int k);

int __inline__ getIndex(int P, int i, int j, int k, int *Index);

void setIndex_arr(int P, int N, int *indices, int indicesSize, int *ii, int iiSize, int *jj, int jjSize, int *kk, int kkSize);

void getCoeff(REAL *__restrict__ a, int aSize, REAL *__restrict__ ax, int axSize, REAL *__restrict__ ay, int aySize, REAL *__restrict__ az, int azSize, REAL dx, REAL dy, REAL dz, int P, REAL kappa, int *__restrict__ Index, int indSize);

void M2P_c(REAL *p, int pSize, REAL *px, int pxSize, REAL *py, int pySize, REAL *pz, int pzSize, REAL *mult, int multSize, REAL *x, int xSize, REAL *y, int ySize, REAL *z, int zSize, REAL *xc, int xcSize, REAL *yc, int ycSize, REAL *zc, int zcSize, int *offMlt, int offMltSize, int *offTar, int offTarSize, int *sizeTar, int sizeTarSize, int P, REAL kappa, int Nm, int *Index, int indSize);



double get_time (void)
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return (double)(tv.tv_sec+1e-6*tv.tv_usec);
};

int setIndex(int P, int i, int j, int k)
{

    int I=0, ii, jj;
    for (ii=0; ii<i; ii++)
    {
        for (jj=1; jj<P+2-ii; jj++)
        {
            I+=jj;
        }
    }
    for (jj=P+2-j; jj<P+2; jj++)
    {
        I+=jj-i;
    }
    I+=k;

    return I;
};

int __inline__ getIndex(int P, int i, int j, int k, int *Index)
{
    return Index[i*(P+1)*(P+1)+j*(P+1)+k];
};

void setIndex_arr(int P, int N, int *indices, int indicesSize, int *ii, int iiSize, int *jj, int jjSize, int *kk, int kkSize)
{
    for (int iter=0; iter<N; iter++)
        indices[iter] = setIndex(P, ii[iter], jj[iter], kk[iter]);

};

void getCoeff(REAL *__restrict__ a, int aSize, REAL *__restrict__ ax, int axSize, REAL *__restrict__ ay, int aySize, REAL *__restrict__ az, int azSize, REAL dx, REAL dy, REAL dz, int P, REAL kappa, int *__restrict__ Index, int indSize)
{
    
    REAL b[aSize];
    REAL R2 = dx*dx+dy*dy+dz*dz;
    REAL R  = sqrt(R2);
    int i,j,k,I,Im1x,Im2x,Im1y,Im2y,Im1z,Im2z;
    REAL C,C1,C2,Cb, R2_1;

    R2_1 = 1/R2;

    // First coefficient
    b[0] = exp(-kappa*R);
    a[0] = b[0]/R;

    // Two indices = 0
    I = getIndex(P,1,0,0,Index);
    b[I]   = -kappa * (dx*a[0]); // 1,0,0
    b[P+1] = -kappa * (dy*a[0]); // 0,1,0
    b[1]   = -kappa * (dz*a[0]); // 0,0,1

    a[I]   = -R2_1*dx*(kappa*b[0]+a[0]);
    a[P+1] = -R2_1*dy*(kappa*b[0]+a[0]);
    a[1]   = -R2_1*dz*(kappa*b[0]+a[0]);
    
    ax[0]  = a[I]; 
    ay[0]  = a[P+1]; 
    az[0]  = a[1]; 

    for (i=2; i<P+1; i++)
    {
        Cb   = -kappa/i;
        C    = R2_1/i;
        I    = getIndex(P,i,0,0,Index);
        Im1x = getIndex(P,i-1,0,0,Index);
        Im2x = getIndex(P,i-2,0,0,Index);
        b[I] = Cb * (dx*a[Im1x] + a[Im2x]);
        a[I] = C * ( -kappa*(dx*b[Im1x] + b[Im2x]) -(2*i-1)*dx*a[Im1x] - (i-1)*a[Im2x] );
        ax[Im1x] = a[I]*i;

        I    = getIndex(P,0,i,0,Index);
        Im1y = I-(P+2-i);
        Im2y = Im1y-(P+2-i+1);
        b[I] = Cb * (dy*a[Im1y] + a[Im2y]);
        a[I] = C * ( -kappa*(dy*b[Im1y] + b[Im2y]) -(2*i-1)*dy*a[Im1y] - (i-1)*a[Im2y] );
        ay[Im1y] = a[I]*i;

        I   = i;
        Im1z = I-1;
        Im2z = I-2;
        b[I] = Cb * (dz*a[Im1z] + a[Im2z]);
        a[I] = C * ( -kappa*(dz*b[Im1z] + b[Im2z]) -(2*i-1)*dz*a[Im1z] - (i-1)*a[Im2z] );
        az[Im1z] = a[I]*i;
    }

    // One index = 0, one = 1 other >=1

    Cb   = -kappa/2;
    C    = R2_1/2.;
    I    = getIndex(P,1,1,0,Index);
    Im1x = P+1;
    Im1y = I-P;
    b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y]);
    a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]) - 3*(dx*a[Im1x]+dy*a[Im1y]) );
    ax[Im1x] = a[I];
    ay[Im1y] = a[I];
    I    = getIndex(P,1,0,1,Index);
    Im1x = 1;
    Im1z = I-1;
    b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z]);
    a[I] = C * ( -kappa*(dx*b[Im1x]+dz*b[Im1z]) - 3*(dx*a[Im1x]+dz*a[Im1z]) );
    ax[Im1x] = a[I];
    az[Im1z] = a[I];
    I    = getIndex(P,0,1,1,Index);
    Im1y = I-(P+1);
    Im1z = I-1;
    b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z]);
    a[I] = C * ( -kappa*(dy*b[Im1y]+dz*b[Im1z]) - 3*(dy*a[Im1y]+dz*a[Im1z]) );
    ay[Im1y] = a[I];
    az[Im1z] = a[I];

    for (i=2; i<P; i++)
    {
        Cb   = -kappa/(i+1);
        C    = R2_1/(1+i);
        C1   = 1+2*i;
        I    = getIndex(P,1,i,0,Index);
        Im1x = getIndex(P,0,i,0,Index);
        Im1y = I-(P+1-i);
        Im2y = Im1y-(P+2-i);
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + a[Im2y]);
        a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+b[Im2y]) - C1*(dx*a[Im1x]+dy*a[Im1y]) - (1+i-1)*(a[Im2y]) );
        ax[Im1x] = a[I];
        ay[Im1y] = a[I]*i;

        I    = getIndex(P,1,0,i,Index);
        Im1x = getIndex(P,0,0,i,Index);
        Im1z = I-1;
        Im2z = I-2;
        b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z] + a[Im2z]);
        a[I] = C * ( -kappa*(dx*b[Im1x]+dz*b[Im1z]+b[Im2z]) - C1*(dx*a[Im1x]+dz*a[Im1z]) - (1+i-1)*(a[Im2z]) );
        ax[Im1x] = a[I];
        az[Im1z] = a[I]*i;
        
        I    = getIndex(P,0,1,i,Index);
        Im1y = I-(P+1);
        Im1z = I-1;
        Im2z = I-2;
        b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z] + a[Im2z]);
        a[I] = C * ( -kappa*(dy*b[Im1y]+dz*b[Im1z]+b[Im2z]) - C1*(dy*a[Im1y]+dz*a[Im1z]) - (1+i-1)*(a[Im2z]) );
        ay[Im1y] = a[I];
        az[Im1z] = a[I]*i;
        
        I    = getIndex(P,i,1,0,Index);
        Im1y = I-(P+1-i);
        Im1x = getIndex(P,i-1,1,0,Index);
        Im2x = getIndex(P,i-2,1,0,Index);
        b[I] = Cb * (dy*a[Im1y] + dx*a[Im1x] + a[Im2x]);
        a[I] = C * ( -kappa*(dy*b[Im1y]+dx*b[Im1x]+b[Im2x]) - C1*(dy*a[Im1y]+dx*a[Im1x]) - (1+i-1)*(a[Im2x]) );
        ax[Im1x] = a[I]*i;
        ay[Im1y] = a[I];
        
        I    = getIndex(P,i,0,1,Index);
        Im1z = I-1;
        Im1x = getIndex(P,i-1,0,1,Index);
        Im2x = getIndex(P,i-2,0,1,Index);
        b[I] = Cb * (dz*a[Im1z] + dx*a[Im1x] + a[Im2x]);
        a[I] = C * ( -kappa*(dz*b[Im1z]+dx*b[Im1x]+b[Im2x]) - C1*(dz*a[Im1z]+dx*a[Im1x]) - (1+i-1)*(a[Im2x]) );
        ax[Im1x] = a[I]*i;
        az[Im1z] = a[I];

        I    = getIndex(P,0,i,1,Index);
        Im1z = I-1;
        Im1y = I-(P+2-i);
        Im2y = Im1y-(P+3-i);
        b[I] = Cb * (dz*a[Im1z] + dy*a[Im1y] + a[Im2y]);
        a[I] = C * ( -kappa*(dz*b[Im1z]+dy*b[Im1y]+b[Im2y]) - C1*(dz*a[Im1z]+dy*a[Im1y]) - (1+i-1)*(a[Im2y]) );
        ay[Im1y] = a[I]*i;
        az[Im1z] = a[I];
    }

    // One index 0, others >=2
    for (i=2; i<P+1; i++)
    {
        for (j=2; j<P+1-i; j++)
        {
            Cb   = -kappa/(i+j); 
            C    = R2_1/(i+j);
            C1   = 2*(i+j)-1;
            I    = getIndex(P,i,j,0,Index);
            Im1x = getIndex(P,i-1,j,0,Index);
            Im2x = getIndex(P,i-2,j,0,Index);
            Im1y = I-(P+2-j-i);
            Im2y = Im1y-(P+3-j-i);
            b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + a[Im2x] + a[Im2y]);
            a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+b[Im2x]+b[Im2y]) - C1*(dx*a[Im1x]+dy*a[Im1y]) -(i+j-1)*(a[Im2x]+a[Im2y]) );
            ax[Im1x] = a[I]*i;
            ay[Im1y] = a[I]*j;

            I    = getIndex(P,i,0,j,Index);
            Im1x = getIndex(P,i-1,0,j,Index);
            Im2x = getIndex(P,i-2,0,j,Index);
            Im1z = I-1;
            Im2z = I-2;
            b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z] + a[Im2x] + a[Im2z]);
            a[I] = C * ( -kappa*(dx*b[Im1x]+dz*b[Im1z]+b[Im2x]+b[Im2z]) - C1*(dx*a[Im1x]+dz*a[Im1z]) -(i+j-1)*(a[Im2x]+a[Im2z]) );
            ax[Im1x] = a[I]*i;
            az[Im1z] = a[I]*j;

            I    = getIndex(P,0,i,j,Index);
            Im1y = I-(P+2-i);
            Im2y = Im1y-(P+3-i);
            Im1z = I-1;
            Im2z = I-2; 
            b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z] + a[Im2y] + a[Im2z]);
            a[I] = C * ( -kappa*(dy*b[Im1y]+dz*b[Im1z]+b[Im2y]+b[Im2z]) - C1*(dy*a[Im1y]+dz*a[Im1z]) -(i+j-1)*(a[Im2y]+a[Im2z]) );
            ay[Im1y] = a[I]*i;
            az[Im1z] = a[I]*j;
        }
    }

    if (P>2)
    {
        // Two index = 1, other>=1
        C    = R2_1/3;
        Cb   = -kappa/3;
        I    = getIndex(P,1,1,1,Index);
        Im1x = getIndex(P,0,1,1,Index);
        Im1y = getIndex(P,1,0,1,Index);
        Im1y = I-(P);
        Im1z = I-1;
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z]);
        a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]) - 5*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) );
        ax[Im1x] = a[I];
        ay[Im1y] = a[I];
        az[Im1z] = a[I];
        for (i=2; i<P-1; i++)
        {
            Cb   = -kappa/(2+i);
            C    = R2_1/(i+2);
            C1   = 2*i+3;
            I    = getIndex(P,i,1,1,Index);
            Im1x = getIndex(P,i-1,1,1,Index);
            Im1y = I-(P+1-i);
            Im1z = I-1;
            Im2x = getIndex(P,i-2,1,1,Index);
            b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x]);
            a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]) - C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2x]) );
            ax[Im1x] = a[I]*i;
            ay[Im1y] = a[I];
            az[Im1z] = a[I];

            I    = getIndex(P,1,i,1,Index);
            Im1x = getIndex(P,0,i,1,Index);
            Im1y = I-(P+1-i);
            Im2y = Im1y-(P+2-i);
            Im1z = I-1 ;
            b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2y]);
            a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2y]) - C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2y]) );
            ax[Im1x] = a[I];
            ay[Im1y] = a[I]*i;
            az[Im1z] = a[I];
            
            I    = getIndex(P,1,1,i,Index);
            Im1x = getIndex(P,0,1,i,Index);
            Im1y = I-(P);
            Im1z = I-1;
            Im2z = I-2;
            b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2z]);
            a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2z]) - C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2z]) );
            ax[Im1x] = a[I];
            ay[Im1y] = a[I];
            az[Im1z] = a[I]*i;
        }
    }

    // One index = 1, others >=2
    if (P>4)
    {
        for (i=2; i<P-2; i++)
        {
            for (j=2; j<P-i; j++)
            {
                Cb = -kappa/(1+i+j);
                C  =  R2_1/(1+i+j);
                C1 = -(2.*(i+j)+1); 
                C2 = (i+j);
                I    = getIndex(P,1,i,j,Index);
                Im1x = getIndex(P,0,i,j,Index);
                Im1y = I-(P+1-i);
                Im2y = Im1y-(P+2-i);
                Im1z = I-1;
                Im2z = I-2; 
                b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2y] + a[Im2z]);
                a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2y]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2y]+a[Im2z]) );
                ax[Im1x] = a[I];
                ay[Im1y] = a[I]*i;
                az[Im1z] = a[I]*j;

                I    = getIndex(P,i,1,j,Index);
                Im1x = getIndex(P,i-1,1,j,Index);
                Im1y = I-(P+1-i);
                Im2x = getIndex(P,i-2,1,j,Index);
                Im1z = I-1;
                Im2z = I-2;
                b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2z]);
                a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2z]) );
                ax[Im1x] = a[I]*i;
                ay[Im1y] = a[I];
                az[Im1z] = a[I]*j;
                
                I    = getIndex(P,i,j,1,Index);
                Im1x = getIndex(P,i-1,j,1,Index);
                Im2x = getIndex(P,i-2,j,1,Index);
                Im1y = I-(P+2-i-j);
                Im2y = Im1y-(P+3-i-j); 
                Im1z = I-1;
                b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2y]);
                a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2y]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2y]) );
                ax[Im1x] = a[I]*i;
                ay[Im1y] = a[I]*j;
                az[Im1z] = a[I];
            }
        }
    }

    // All indices >= 2
    if (P>5)
    {
        for (i=2;i<P-3;i++)
        {
            for (j=2;j<P-1-i;j++)
            {
                for (k=2;k<P+1-i-j;k++)
                {
                    Cb = -kappa/(i+j+k);
                    C  = R2_1/(i+j+k);
                    C1 = -(2.*(i+j+k)-1); 
                    C2 = i+j+k-1.;
                    I    = getIndex(P,i,j,k,Index);
                    Im1x = getIndex(P,i-1,j,k,Index);
                    Im2x = getIndex(P,i-2,j,k,Index);
                    Im1y = I-(P+2-i-j);
                    Im2y = Im1y-(P+3-i-j); 
                    Im1z = I-1; 
                    Im2z = I-2; 
                    b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2y] + a[Im2z]);
                    a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2y]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2y]+a[Im2z]) );
                    ax[Im1x] = a[I]*i;
                    ay[Im1y] = a[I]*j;
                    az[Im1z] = a[I]*k;
                }
            }
        }
    }
};

void M2P_c(REAL *p, int pSize, REAL *px, int pxSize, REAL *py, int pySize, REAL *pz, int pzSize, 
           REAL *mult, int multSize, REAL *x, int xSize, REAL *y, int ySize, REAL *z, int zSize, 
           REAL *xc, int xcSize, REAL *yc, int ycSize, REAL *zc, int zcSize, int *offMlt, int offMltSize, 
           int *offTar, int offTarSize, int *sizeTar, int sizeTarSize, int P, REAL kappa, int Nm, int *Index, int indSize)
{
    REAL a_aux[Nm];
    REAL ax_aux[Nm];
    REAL ay_aux[Nm];
    REAL az_aux[Nm];

    int CI_start, CI_end, CJ_start, CJ_end;
    REAL dx, dy, dz, sum_p, sum_px, sum_py, sum_pz;

    for (int ii=0; ii<Nm; ii++)
    {
        a_aux[ii] = 0.;
        ax_aux[ii] = 0.;
        ay_aux[ii] = 0.;
        az_aux[ii] = 0.;
    }

    for (int tarTwg=0; tarTwg<offTarSize; tarTwg++)
    {
        CI_start = offTar[tarTwg];
        CI_end = CI_start + sizeTar[tarTwg];
        CJ_start = offMlt[tarTwg];
        CJ_end = offMlt[tarTwg+1];

        for (int i=CI_start; i<CI_end; i++)
        {
            sum_p = 0.;
            sum_px = 0.;
            sum_py = 0.;
            sum_pz = 0.;

            for (int j=CJ_start; j<CJ_end; j++)
            {
                dx = x[i] - xc[j];
                dy = y[i] - yc[j];
                dz = z[i] - zc[j];

                getCoeff(a_aux, Nm, ax_aux, Nm, ay_aux, Nm, az_aux, Nm, dx, dy, dz, P, kappa, Index, indSize);

                for (int jj=0; jj<Nm; jj++)
                {
                    sum_p += a_aux[jj]*mult[Nm*j+jj];
                    sum_px += ax_aux[jj]*mult[Nm*j+jj];
                    sum_py += ay_aux[jj]*mult[Nm*j+jj];
                    sum_pz += az_aux[jj]*mult[Nm*j+jj];
                }
            }

            p[i]  += sum_p;
            px[i] += sum_px;
            py[i] += sum_py;
            pz[i] += sum_pz;
        }
    }
}
