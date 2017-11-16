from pycuda.compiler import SourceModule

def kernels(BSZ, Nm, P):
    
    mod = SourceModule( """

    #define REAL double
    #define BSZ %(blocksize)d
    #define Nm  %(Nmult)d
    #define P   %(Ptree)d

    /*
    __device__ int getIndex(int i, int j, int k)
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
    }
    */
    __device__ int getIndex(int i, int j, int k, int *Index)
    {   
        return Index[(P+1)*(P+1)*i + (P+1)*j + k]; 
    }


    __device__ void getCoeff(REAL *a, REAL *ax, REAL *ay, REAL *az, 
                            REAL dx, REAL dy, REAL dz, REAL kappa, int *Index)
    {
        REAL b[Nm];
        REAL R = sqrt(dx*dx+dy*dy+dz*dz);
        REAL R2 = R*R;
        
        int i,j,k,I,Im1x,Im2x,Im1y,Im2y,Im1z,Im2z;
        REAL C,C1,C2,Cb;

        // First coefficient
        
        b[0] = exp(-kappa*R);
        a[0] = b[0]/R;

        // Two indices = 0
        I = getIndex(1,0,0, Index);
        b[I]   = -kappa * (dx*a[0]); // 1,0,0
        b[P+1] = -kappa * (dy*a[0]); // 0,1,0
        b[1]   = -kappa * (dz*a[0]); // 0,0,1

        a[I]   = -1/R2*(kappa*dx*b[0]+dx*a[0]);
        a[P+1] = -1/R2*(kappa*dy*b[0]+dy*a[0]);
        a[1]   = -1/R2*(kappa*dz*b[0]+dz*a[0]);
     
        ax[0]  = a[I]; 
        ay[0]  = a[P+1]; 
        az[0]  = a[1]; 
      
        for (i=2; i<P+1; i++)
        {   
            Cb   = -kappa/i;
            C    = 1./(i*R2);
            I    = getIndex(i,0,0, Index);
            Im1x = getIndex(i-1,0,0, Index);
            Im2x = getIndex(i-2,0,0, Index);
            b[I] = Cb * (dx*a[Im1x] + a[Im2x]);
            a[I] = C * ( -kappa*(dx*b[Im1x] + b[Im2x]) -(2*i-1)*dx*a[Im1x] - (i-1)*a[Im2x] );
            ax[Im1x] = a[I]*i;

            I    = getIndex(0,i,0, Index);
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
        I    = getIndex(1,1,0, Index);
        Im1x = P+1;
        Im1y = I-(P+2-1-1);
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y]);
        a[I] = 1./(2*R2) * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]) -(2*2-1)*(dx*a[Im1x]+dy*a[Im1y]) );
        ax[Im1x] = a[I];
        ay[Im1y] = a[I];
        I    = getIndex(1,0,1, Index);
        Im1x = 1;
        Im1z = I-1;
        b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z]);
        a[I] = 1./(2*R2) * ( -kappa*(dx*b[Im1x]+dz*b[Im1z]) -(2*2-1)*(dx*a[Im1x]+dz*a[Im1z]) );
        ax[Im1x] = a[I];
        az[Im1z] = a[I];
        I    = getIndex(0,1,1, Index);
        Im1y = I-(P+2-1);
        Im1z = I-1;
        b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z]);
        a[I] = 1./(2*R2) * ( -kappa*(dy*b[Im1y]+dz*b[Im1z]) -(2*2-1)*(dy*a[Im1y]+dz*a[Im1z]) );
        ay[Im1y] = a[I];
        az[Im1z] = a[I];

        for (i=2; i<P; i++)
        {
            Cb   = -kappa/(i+1);
            C    = 1./((1+i)*R2);
            I    = getIndex(1,i,0, Index);
            Im1x = getIndex(0,i,0, Index);
            Im1y = I-(P+2-i-1);
            Im2y = Im1y-(P+2-i);
            b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + a[Im2y]);
            a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+b[Im2y]) -(2*(1+i)-1)*(dx*a[Im1x]+dy*a[Im1y]) - (1+i-1)*(a[Im2y]) );
            ax[Im1x] = a[I];
            ay[Im1y] = a[I]*i;

            I    = getIndex(1,0,i, Index);
            Im1x = getIndex(0,0,i, Index);
            Im1z = I-1;
            Im2z = I-2;
            b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z] + a[Im2z]);
            a[I] = C * ( -kappa*(dx*b[Im1x]+dz*b[Im1z]+b[Im2z]) -(2*(1+i)-1)*(dx*a[Im1x]+dz*a[Im1z]) - (1+i-1)*(a[Im2z]) );
            ax[Im1x] = a[I];
            az[Im1z] = a[I]*i;

            I    = getIndex(0,1,i, Index);
            Im1y = I-(P+2-1);
            Im1z = I-1;
            Im2z = I-2;
            b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z] + a[Im2z]);
            a[I] = C * ( -kappa*(dy*b[Im1y]+dz*b[Im1z]+b[Im2z]) -(2*(1+i)-1)*(dy*a[Im1y]+dz*a[Im1z]) - (1+i-1)*(a[Im2z]) );
            ay[Im1y] = a[I];
            az[Im1z] = a[I]*i;

            I    = getIndex(i,1,0, Index);
            Im1y = I-(P+2-1-i);
            Im1x = getIndex(i-1,1,0, Index);
            Im2x = getIndex(i-2,1,0, Index);
            b[I] = Cb * (dy*a[Im1y] + dx*a[Im1x] + a[Im2x]);
            a[I] = C * ( -kappa*(dy*b[Im1y]+dx*b[Im1x]+b[Im2x]) -(2*(1+i)-1)*(dy*a[Im1y]+dx*a[Im1x]) - (1+i-1)*(a[Im2x]) );
            ax[Im1x] = a[I]*i;
            ay[Im1y] = a[I];

            I    = getIndex(i,0,1, Index);
            Im1z = I-1;
            Im1x = getIndex(i-1,0,1, Index);
            Im2x = getIndex(i-2,0,1, Index);
            b[I] = Cb * (dz*a[Im1z] + dx*a[Im1x] + a[Im2x]);
            a[I] = C * ( -kappa*(dz*b[Im1z]+dx*b[Im1x]+b[Im2x]) -(2*(1+i)-1)*(dz*a[Im1z]+dx*a[Im1x]) - (1+i-1)*(a[Im2x]) );
            ax[Im1x] = a[I]*i;
            az[Im1z] = a[I];

            I    = getIndex(0,i,1, Index);
            Im1z = I-1;
            Im1y = I-(P+2-i);
            Im2y = Im1y-(P+2-i+1);
            b[I] = Cb * (dz*a[Im1z] + dy*a[Im1y] + a[Im2y]);
            a[I] = C * ( -kappa*(dz*b[Im1z]+dy*b[Im1y]+b[Im2y]) -(2*(1+i)-1)*(dz*a[Im1z]+dy*a[Im1y]) - (1+i-1)*(a[Im2y]) );
            ay[Im1y] = a[I]*i;
            az[Im1z] = a[I];
        }

        // One index 0, others >=2
        for (i=2; i<P+1; i++)
        {
            for (j=2; j<P+1-i; j++)
            {
                Cb   = -kappa/(i+j);
                C    = 1./((i+j)*R2);
                I    = getIndex(i,j,0, Index);
                Im1x = getIndex(i-1,j,0, Index);
                Im2x = getIndex(i-2,j,0, Index);
                Im1y = I-(P+2-j-i);
                Im2y = Im1y-(P+3-j-i);
                b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + a[Im2x] + a[Im2y]);
                a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+b[Im2x]+b[Im2y]) -(2*(i+j)-1)*(dx*a[Im1x]+dy*a[Im1y]) -(i+j-1)*(a[Im2x]+a[Im2y]) );
                ax[Im1x] = a[I]*i;
                ay[Im1y] = a[I]*j;

                I    = getIndex(i,0,j, Index);
                Im1x = getIndex(i-1,0,j, Index);
                Im2x = getIndex(i-2,0,j, Index);
                Im1z = I-1;
                Im2z = I-2;
                b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z] + a[Im2x] + a[Im2z]);
                a[I] = C * ( -kappa*(dx*b[Im1x]+dz*b[Im1z]+b[Im2x]+b[Im2z]) -(2*(i+j)-1)*(dx*a[Im1x]+dz*a[Im1z]) -(i+j-1)*(a[Im2x]+a[Im2z]) );
                ax[Im1x] = a[I]*i;
                az[Im1z] = a[I]*j;

                I    = getIndex(0,i,j, Index);
                Im1y = I-(P+2-i);
                Im2y = Im1y-(P+3-i);
                Im1z = I-1;
                Im2z = I-2;
                b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z] + a[Im2y] + a[Im2z]);
                a[I] = C * ( -kappa*(dy*b[Im1y]+dz*b[Im1z]+b[Im2y]+b[Im2z]) -(2*(i+j)-1)*(dy*a[Im1y]+dz*a[Im1z]) -(i+j-1)*(a[Im2y]+a[Im2z]) );
                ay[Im1y] = a[I]*i;
                az[Im1z] = a[I]*j;
            }
        }

        if (P>2)
        {
            // Two index = 1, other>=1
            Cb   = -kappa/3;
            I    = getIndex(1,1,1, Index);
            Im1x = getIndex(0,1,1, Index);
            Im1y = getIndex(1,0,1, Index);
            Im1y = I-(P);
            Im1z = I-1;
            b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z]);
            a[I] = 1/(3*R2) * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]) -5*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) );
            ax[Im1x] = a[I];
            ay[Im1y] = a[I];
            az[Im1z] = a[I];

            for (i=2; i<P-1; i++)
            {
                Cb   = -kappa/(2+i);
                C    = 1./((i+2)*R2);
                I    = getIndex(i,1,1, Index);
                Im1x = getIndex(i-1,1,1, Index);
                Im1y = I-(P+2-i-1);
                Im1z = I-1;
                Im2x = getIndex(i-2,1,1, Index);
                b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x]);
                a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]) -(2*(i+2)-1)*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2x]) );
                ax[Im1x] = a[I]*i;
                ay[Im1y] = a[I];
                az[Im1z] = a[I];

                I    = getIndex(1,i,1, Index);
                Im1x = getIndex(0,i,1, Index);
                Im1y = I-(P+2-i-1);
                Im2y = Im1y-(P+3-i-1);
                Im1z = I-1 ;
                b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2y]);
                a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2y]) -(2*(i+2)-1)*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2y]) );
                ax[Im1x] = a[I];
                ay[Im1y] = a[I]*i;
                az[Im1z] = a[I];

                I    = getIndex(1,1,i, Index);
                Im1x = getIndex(0,1,i, Index);
                Im1y = I-(P);
                Im1z = I-1;
                Im2z = I-2;
                b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2z]);
                a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2z]) -(2*(i+2)-1)*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2z]) );
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
                    C  = 1./((1+i+j)*R2);
                    C1 = -(2.*(1+i+j)-1);
                    C2 = (i+j);
                    I    = getIndex(1,i,j, Index);
                    Im1x = getIndex(0,i,j, Index);
                    Im1y = I-(P+2-1-i);
                    Im2y = Im1y-(P+3-1-i);
                    Im1z = I-1;
                    Im2z = I-2;
                    b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2y] + a[Im2z]);
                    a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2y]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2y]+a[Im2z]) );
                    ax[Im1x] = a[I];
                    ay[Im1y] = a[I]*i;
                    az[Im1z] = a[I]*j;

                    I    = getIndex(i,1,j, Index);
                    Im1x = getIndex(i-1,1,j, Index);
                    Im1y = I-(P+2-i-1);
                    Im2x = getIndex(i-2,1,j, Index);
                    Im1z = I-1;
                    Im2z = I-2;
                    b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2z]);
                    a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2z]) );
                    ax[Im1x] = a[I]*i;
                    ay[Im1y] = a[I];
                    az[Im1z] = a[I]*j;

                    I    = getIndex(i,j,1, Index);
                    Im1x = getIndex(i-1,j,1, Index);
                    Im2x = getIndex(i-2,j,1, Index);
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
                        C  = 1./((i+j+k)*R2);
                        C1 = -(2.*(i+j+k)-1);
                        C2 = i+j+k-1.;
                        I    = getIndex(i,j,k, Index);
                        Im1x = getIndex(i-1,j,k, Index);
                        Im2x = getIndex(i-2,j,k, Index);
                        Im1y = I-(P+2-i-j);
                        Im2y = Im1y-(P+3-i-j);
                        Im1z = I-1;
                        Im2z = I-2;
                        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2y] + a[Im2z]);
                        a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2y]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2y]+a[Im2z]) );
                        ax[Im1x] = a[I]*i;
                        ay[Im1y] = a[I]*j;
                        az[Im1z] = a[I]*k;

                        ax[I] = 0.;
                        ay[I] = 0.;
                        az[I] = 0.;
                    }
                }
            }
        }

    }

    __device__ void multipole(REAL &s, REAL &sx, REAL &sy, REAL &sz, 
                            REAL *mp, int CJ_start, int jblock, 
                            REAL *coeff, REAL *coeff_x, REAL *coeff_y, REAL *coeff_z, int j)
    {

        for (int i_coeff=0; i_coeff<Nm; i_coeff++)
        {
            s  += mp[(CJ_start+j)*Nm + jblock*BSZ*Nm + i_coeff] * coeff[i_coeff];
            sx += mp[(CJ_start+j)*Nm + jblock*BSZ*Nm + i_coeff] * coeff_x[i_coeff];
            sy += mp[(CJ_start+j)*Nm + jblock*BSZ*Nm + i_coeff] * coeff_y[i_coeff];
            sz += mp[(CJ_start+j)*Nm + jblock*BSZ*Nm + i_coeff] * coeff_z[i_coeff];
        }

    }

    __global__ void M2P(int *sizeTarDev, int *offsetMltDev,
                        REAL *xtDev, REAL *ytDev, REAL *ztDev,
                        REAL *xcDev, REAL *ycDev, REAL *zcDev, REAL *mpDev,
                        REAL *p, REAL *px, REAL *py, REAL *pz,
                        REAL kappa, int BpT, int *IndexDev, int NCRIT)
    {
        int I = threadIdx.x + blockIdx.x*NCRIT;
        int CJ_start = offsetMltDev[blockIdx.x];
        int Nmlt     = offsetMltDev[blockIdx.x+1] - CJ_start;

        REAL xi, yi, zi;
        REAL dx, dy, dz, s = 0., sx = 0., sy = 0., sz = 0.;
        REAL coeff[Nm], coeff_x[Nm], coeff_y[Nm], coeff_z[Nm];

        __shared__ REAL xc_sh[BSZ],
                        yc_sh[BSZ],
                        zc_sh[BSZ];
        __shared__ int Index_sh[(P+1)*(P+1)*(P+1)]; 
        int i;

        for (int ind=0; ind<((P+1)*(P+1)*(P+1)-1)/BSZ; ind++)
        {
            Index_sh[ind*BSZ + threadIdx.x] = IndexDev[ind*BSZ + threadIdx.x];    
        }

        int ind = ((P+1)*(P+1)*(P+1)-1)/BSZ;
        if (threadIdx.x<(P+1)*(P+1)*(P+1)-BSZ*ind)
        {
            Index_sh[ind*BSZ + threadIdx.x] = IndexDev[ind*BSZ + threadIdx.x];
        }

        __syncthreads();

        for (int mult=0; mult<Nm; mult++)
        {
            coeff_x[mult] = 0.;
            coeff_y[mult] = 0.;
            coeff_z[mult] = 0.;
        }

        for (int iblock=0; iblock<BpT; iblock++)
        {
            i  = I + iblock*BSZ;
            xi = xtDev[i];
            yi = ytDev[i];
            zi = ztDev[i];

            for(int jblock=0; jblock<(Nmlt-1)/BSZ; jblock++)
            {
                __syncthreads();
                xc_sh[threadIdx.x] = xcDev[CJ_start + jblock*BSZ + threadIdx.x];
                yc_sh[threadIdx.x] = ycDev[CJ_start + jblock*BSZ + threadIdx.x];
                zc_sh[threadIdx.x] = zcDev[CJ_start + jblock*BSZ + threadIdx.x];
                __syncthreads();

                if (threadIdx.x+iblock*BSZ<sizeTarDev[blockIdx.x])
                {
                    for (int j=0; j<BSZ; j++)
                    {
                        dx = xi - xc_sh[j];
                        dy = yi - yc_sh[j];
                        dz = zi - zc_sh[j];
                        getCoeff(coeff, coeff_x, coeff_y, coeff_z, dx, dy, dz, kappa, Index_sh);
                        multipole(s , sx, sy, sz, mpDev, CJ_start, jblock, coeff, coeff_x, coeff_y, coeff_z, j);
                    }
                }
            } 

            __syncthreads();
            int jblock = (Nmlt-1)/BSZ;
            xc_sh[threadIdx.x] = xcDev[CJ_start + jblock*BSZ + threadIdx.x];
            yc_sh[threadIdx.x] = ycDev[CJ_start + jblock*BSZ + threadIdx.x];
            zc_sh[threadIdx.x] = zcDev[CJ_start + jblock*BSZ + threadIdx.x];
            __syncthreads();
            
            if (threadIdx.x+iblock*BSZ<sizeTarDev[blockIdx.x])
            {
                for (int j=0; j<Nmlt-(jblock*BSZ); j++)
                {
                    dx = xi - xc_sh[j];
                    dy = yi - yc_sh[j];
                    dz = zi - zc_sh[j];
                    getCoeff(coeff, coeff_x, coeff_y, coeff_z, dx, dy, dz, kappa, Index_sh);
                    multipole(s , sx, sy, sz, mpDev, CJ_start, jblock, coeff, coeff_x, coeff_y, coeff_z, j);
                }
            }

            if (threadIdx.x+iblock*BSZ<sizeTarDev[blockIdx.x])
            {
                p[i] += s;
                px[i] += sx;
                py[i] += sy;
                pz[i] += sz;
            }
        }
        
    }


    __global__ void P2P(int *offsetSrcDev, int *offsetIntDev, int *intPtrDev, int *sizeTarDev,
                        REAL *xsDev, REAL *ysDev, REAL *zsDev, REAL *mDev,
                        REAL *xtDev, REAL *ytDev, REAL *ztDev, REAL *p,
                        REAL *px, REAL *py, REAL *pz, REAL kappa, REAL eps, int BpT, int NCRIT)
    {
        int I = threadIdx.x + blockIdx.x*NCRIT;
        int CJt_start = offsetIntDev[blockIdx.x];
        int Ntsrc     = offsetIntDev[blockIdx.x+1] - CJt_start;
        
        REAL xi, yi, zi;

        REAL dx, dy, dz, r, s, sx, sy, sz;

        __shared__ REAL xj_sh[BSZ],
                        yj_sh[BSZ],
                        zj_sh[BSZ],
                        m_sh[BSZ];
        int i, CJ_start, Nsrc, twig_ptr;

        for (int iblock=0; iblock<BpT; iblock++)
        {
            i  = I + iblock*BSZ;
            xi = xtDev[i];
            yi = ytDev[i];
            zi = ztDev[i];
            s = 0., sx= 0., sy=0., sz=0.;

            for (int twig=0; twig<Ntsrc; twig++)
            {
                twig_ptr = intPtrDev[CJt_start + twig];
                CJ_start = offsetSrcDev[twig_ptr];
                Nsrc = offsetSrcDev[twig_ptr+1] - CJ_start;

                for(int jblock=0; jblock<(Nsrc-1)/BSZ; jblock++)
                {
                    __syncthreads();
                    xj_sh[threadIdx.x] = xsDev[CJ_start + jblock*BSZ + threadIdx.x];
                    yj_sh[threadIdx.x] = ysDev[CJ_start + jblock*BSZ + threadIdx.x];
                    zj_sh[threadIdx.x] = zsDev[CJ_start + jblock*BSZ + threadIdx.x];
                    m_sh[threadIdx.x] = mDev[CJ_start + jblock*BSZ + threadIdx.x];
                    __syncthreads();
                    
                    if (threadIdx.x+iblock*BSZ<sizeTarDev[blockIdx.x])
                    {
                        for (int j=0; j<BSZ; j++)
                        {
                            dx = xi - xj_sh[j];
                            dy = yi - yj_sh[j];
                            dz = zi - zj_sh[j];
                            r = sqrt(dx*dx + dy*dy + dz*dz + eps);
                            s  += m_sh[j] * exp(-kappa*r)/r; 
                            sx += -m_sh[j] * dx * exp(-kappa*r)/(r*r)*(kappa+1/r); 
                            sy += -m_sh[j] * dy * exp(-kappa*r)/(r*r)*(kappa+1/r); 
                            sz += -m_sh[j] * dz * exp(-kappa*r)/(r*r)*(kappa+1/r); 
                        }
                    }
                }
               
                __syncthreads();
                int jblock = (Nsrc-1)/BSZ;
                xj_sh[threadIdx.x] = xsDev[CJ_start + jblock*BSZ + threadIdx.x];
                yj_sh[threadIdx.x] = ysDev[CJ_start + jblock*BSZ + threadIdx.x];
                zj_sh[threadIdx.x] = zsDev[CJ_start + jblock*BSZ + threadIdx.x];
                m_sh[threadIdx.x] = mDev[CJ_start + jblock*BSZ + threadIdx.x];
                __syncthreads();

                if (threadIdx.x+iblock*BSZ<sizeTarDev[blockIdx.x])
                {
                    for (int j=0; j<Nsrc-(jblock*BSZ); j++)
                    {
                        dx = xi - xj_sh[j];
                        dy = yi - yj_sh[j];
                        dz = zi - zj_sh[j];
                        r = sqrt(dx*dx + dy*dy + dz*dz + eps);
                        s  += m_sh[j] * exp(-kappa*r)/r; 
                        sx += -m_sh[j] * dx * exp(-kappa*r)/(r*r)*(kappa+1/r); 
                        sy += -m_sh[j] * dy * exp(-kappa*r)/(r*r)*(kappa+1/r); 
                        sz += -m_sh[j] * dz * exp(-kappa*r)/(r*r)*(kappa+1/r); 
                    }
                }
            }

            if (threadIdx.x+iblock*BSZ<sizeTarDev[blockIdx.x])
            {
                p[i]  += s;
                px[i] += sx;
                py[i] += sy;
                pz[i] += sz;
            }
                
        }
    }

    """%{'blocksize':BSZ, 'Nmult':Nm, 'Ptree':P}, nvcc="nvcc")#, options=["-Xptxas=-v"])

    return mod
