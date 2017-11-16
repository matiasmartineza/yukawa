from numpy import *
from numpy import float64 as REAL
from pylab      import *
from tree_utils import Cell, generateTree, getCoeff_cy, getMultipole, upwardSweep, packData, computeIndices, precompute_terms, addSources, addSources2
from scipy.misc import factorial 
from P2Pwrap import direct_c_cy, P2P_c_cy
from M2Pwrap import M2P_c_cy
import time
import sys
import cProfile

# mode can be 'cpu' or 'pycuda'
mode = 'cpu'

# pyCUDA libraries
if mode == 'pycuda':
    import pycuda.autoinit
    import pycuda.gpuarray  as gpuarray
    import pycuda.driver    as cuda
    from cuda_kernels import kernels

P   = 6
eps = 1e-10
Nm  = int32(factorial(P+3)/(6*factorial(P)))

kappa = 1.

NCRIT = 300
theta = 0.5
BSZ = 256

sample_size = 10000 # Sample for error calculation

Nj  = 100000  # Number of sources
N   = 100000  # Number of targets

# Sources
xj = zeros(Nj)  
yj = zeros(Nj) 
zj = zeros(Nj) 
m = zeros(Nj) 
# Targets
xi = zeros(N)   
yi = zeros(N) 
zi = zeros(N) 

sources_points = open('sources_points.txt', 'w')

for j in range(4):
	for i in random_sample(Nj):
		sources_points.write(str(i) + '\n')

for j in range(3):
	for i in random_sample(N):
		sources_points.write(str(i) + '\n')

sources_points.close()


sources_points = open('sources_points.txt', 'r')

for i in range(Nj):
	xj[i] = float(sources_points.readline())
for i in range(Nj):
	yj[i] = float(sources_points.readline())
for i in range(Nj):
	zj[i] = float(sources_points.readline())
for i in range(Nj):
	m[i] = float(sources_points.readline())
for i in range(N):
	xi[i] = float(sources_points.readline())
for i in range(N):
	yi[i] = float(sources_points.readline())
for i in range(N):
	zi[i] = float(sources_points.readline()) 
	

# Sources
#xj = random_sample(Nj)  
#yj = random_sample(Nj)
#zj = random_sample(Nj)
#m = random_sample(Nj)
# Targets
#xi = random_sample(N)   
#yi = random_sample(N)
#zi = random_sample(N)

x0 = array([0.5,0.5,0.5,0.6]) # Position and radius of first cell

print 'Sources: %i'%Nj
print 'Targets: %i'%N
print 'P      : %i'%P
print 'NCRIT  : %i'%NCRIT
print 'Nm     : %i'%Nm
print 'Sample : %i'%sample_size
print 'BSZ    : %i'%BSZ

phi_d  = zeros(N)
phi_xd = zeros(N)
phi_yd = zeros(N)
phi_zd = zeros(N)

# Direct summation
tic = time.time()
direct_c_cy(xi[0:sample_size],yi[0:sample_size],zi[0:sample_size],xj,yj,zj,phi_d[0:sample_size],phi_xd[0:sample_size],phi_yd[0:sample_size],phi_zd,m,kappa,eps)

toc = time.time()
time_direct = toc-tic

start_tree = time.time()

tic = time.time()
Cells = generateTree(xi, yi, zi, NCRIT, Nm, N, x0)
toc = time.time()
time_tree = toc-tic

tic = time.time()
II,JJ,KK,index, index_large = computeIndices(P)

toc = time.time()
time_ind = toc-tic

tic = time.time()
#addSources(xi,yi,zi,Cells,twig)
for j in range(Nj):
    C = 0
    addSources2(xj,yj,zj,j,Cells,C,NCRIT)
toc = time.time()
time_src = toc-tic

tic = time.time()
C = 0
twig = []
twig = getMultipole(Cells, C, xj, yj, zj, m, twig, II, JJ, KK, index, P, NCRIT)
toc = time.time()
time_P2M = toc-tic


rr = zeros(len(Cells))
for i in range(len(Cells)):
    rr[i] = Cells[i].r
Levels = log(Cells[0].r/min(rr))/log(2) + 1

print 'Cells  : %i'%len(Cells)
print 'Twigs  : %i'%len(twig)
print 'Levels : %i'%Levels


tic = time.time()
combII, combJJ, combKK, IImii, JJmjj, KKmkk, index_small, index_ptr = precompute_terms(P, II, JJ, KK)
toc = time.time()
time_pre = toc-tic

tic = time.time()
for C in reversed(range(1,len(Cells))):
    PC = Cells[C].parent
    upwardSweep(Cells,C,PC,P, II, JJ, KK, index, combII, combJJ, combKK, IImii, JJmjj, KKmkk, index_small, index_ptr)
toc = time.time()
time_M2M = toc-tic

# Data packing
tic = time.time()
Nround = len(twig)*NCRIT
Nlist  = Nround*len(twig)

offsetTarHost = zeros(len(twig)  , dtype=int32)
sizeTarHost   = zeros(len(twig)  , dtype=int32)
offsetSrcHost = zeros(len(twig)+1, dtype=int32)
offsetMltHost = zeros(len(twig)+1, dtype=int32)
offsetIntHost = zeros(len(twig)+1, dtype=int32)

tarPtr = zeros(Nround, dtype=int32)
srcPtr = zeros(Nj, dtype=int32)
mltPtr = zeros(len(twig)*len(twig), dtype=int32)
intPtr = zeros(len(twig)*len(twig), dtype=int32)

t = -1
offInt = 0
offMlt = 0
offSrc = 0
for CI in twig:
    t += 1
    offsetIntHost[t] = offInt
    offsetMltHost[t] = offMlt
    offsetSrcHost[t] = offSrc
    CJ = 0
    intPtr, mltPtr, offInt, offMlt = packData(Cells, CJ, CI, intPtr, mltPtr, offInt, offMlt, theta, NCRIT)
    srcPtr[offSrc:offSrc+Cells[CI].nsource] = Cells[CI].source
    offSrc += Cells[CI].nsource
offsetIntHost[-1] = offInt
offsetMltHost[-1] = offMlt
offsetSrcHost[-1] = offSrc

intPtrHost = zeros(offInt, dtype=int32)
for i in range(len(twig)):
    twig_index = where(intPtr[0:offInt]==twig[i])[0]
    intPtrHost[twig_index] = i

xtHost = zeros(Nround)
ytHost = zeros(Nround)
ztHost = zeros(Nround)

twig_n = -1
for CI in twig:
    twig_n += 1
    targets = Cells[CI].target
    ntargets = Cells[CI].ntarget

    offsetTarHost[twig_n] = twig_n*NCRIT
    sizeTarHost[twig_n] = ntargets

    tarPtr[twig_n*NCRIT:twig_n*NCRIT+Cells[CI].ntarget] = targets
    xtHost[twig_n*NCRIT:twig_n*NCRIT+ntargets] = xi[targets]
    ytHost[twig_n*NCRIT:twig_n*NCRIT+ntargets] = yi[targets]
    ztHost[twig_n*NCRIT:twig_n*NCRIT+ntargets] = zi[targets]

xsHost = xj[srcPtr]
ysHost = yj[srcPtr]
zsHost = zj[srcPtr]
mHost = m[srcPtr]

i = -1
xcHost = zeros(offMlt)
ycHost = zeros(offMlt)
zcHost = zeros(offMlt)
mpHost = zeros(offMlt*Nm)


for C in mltPtr[0:offMlt]:
    i += 1
    xcHost[i] = Cells[C].xc
    ycHost[i] = Cells[C].yc
    zcHost[i] = Cells[C].zc
    mpHost[i*Nm:i*Nm+Nm] = Cells[C].M
toc = time.time()
time_pack = toc - tic

tic = time.time()


if mode == 'pycuda':
    tec = cuda.Event()
    tac = cuda.Event()
    tec.record()
if mode == 'cpu':
    tec = time.time()

time_trans = 0

time_coeff = 0
time_P2P   = 0
p = zeros(N)
px = zeros(N)
py = zeros(N)
pz = zeros(N)

if mode == 'pycuda':
    offsetIntDev = gpuarray.to_gpu(offsetIntHost.astype(int32))
    intPtrDev = gpuarray.to_gpu(intPtrHost.astype(int32))

    pDev = gpuarray.zeros(Nround, dtype=REAL)
    pxDev = gpuarray.zeros(Nround, dtype=REAL)
    pyDev = gpuarray.zeros(Nround, dtype=REAL)
    pzDev = gpuarray.zeros(Nround, dtype=REAL)

    offsetSrcDev = gpuarray.to_gpu(offsetSrcHost.astype(int32))
    offsetMltDev = gpuarray.to_gpu(offsetMltHost.astype(int32))
    sizeTarDev   = gpuarray.to_gpu(sizeTarHost.astype(int32))
    xsDev = gpuarray.to_gpu(xsHost.astype(REAL))
    ysDev = gpuarray.to_gpu(ysHost.astype(REAL))
    zsDev = gpuarray.to_gpu(zsHost.astype(REAL))
    mDev  = gpuarray.to_gpu(mHost.astype(REAL))
    xcDev = gpuarray.to_gpu(xcHost.astype(REAL))
    ycDev = gpuarray.to_gpu(ycHost.astype(REAL))
    zcDev = gpuarray.to_gpu(zcHost.astype(REAL))
    mpDev = gpuarray.to_gpu(mpHost.astype(REAL))
    xtDev = gpuarray.to_gpu(xtHost.astype(REAL))
    ytDev = gpuarray.to_gpu(ytHost.astype(REAL))
    ztDev = gpuarray.to_gpu(ztHost.astype(REAL))
    IndexDev = gpuarray.to_gpu(index_large.astype(int32))

BlocksPerTwig = int(ceil(NCRIT/float(BSZ))) 
print 'Block per twig: %i'%BlocksPerTwig

GSZ = int(ceil(float(Nround)/NCRIT))
    
if mode == 'pycuda':
    mod = kernels(BSZ, Nm, P)

if mode == 'pycuda':
    tac.record()
    tac.synchronize()
    time_trans += tec.time_till(tac)*1e-3

if mode == 'cpu':
    tac = time.time()
    time_trans += tac - tec

paux = zeros(Nround, dtype=REAL)
pxaux = zeros(Nround, dtype=REAL)
pyaux = zeros(Nround, dtype=REAL)
pzaux = zeros(Nround, dtype=REAL)

if mode == 'pycuda':
    tec.record()
    M2P_pycuda = mod.get_function("M2P")
    M2P_pycuda(sizeTarDev, offsetMltDev, xtDev, ytDev, ztDev, xcDev, ycDev, zcDev, mpDev, pDev, pxDev, pyDev, pzDev, REAL(kappa), int32(BlocksPerTwig), IndexDev, int32(NCRIT), block=(BSZ,1,1), grid=(GSZ,1)) 
    tac.record()
    tac.synchronize()
    time_M2P = tec.time_till(tac)*1e-3

if mode == 'cpu':
    tec = time.time()
    M2P_c_cy(paux, pxaux, pyaux, pzaux, mpHost, xtHost, ytHost, ztHost, xcHost, ycHost, zcHost, offsetMltHost, offsetTarHost, sizeTarHost, P, kappa, int(Nm), index_large)
    tac = time.time()
    time_M2P = tac - tec

if mode == 'pycuda':
    tec.record()
    P2P_pycuda = mod.get_function("P2P")
    P2P_pycuda(offsetSrcDev, offsetIntDev, intPtrDev, sizeTarDev, xsDev, ysDev, zsDev, mDev, xtDev, ytDev, ztDev, pDev, pxDev, pyDev, pzDev, REAL(kappa), REAL(eps), int32(BlocksPerTwig), int32(NCRIT), block=(BSZ,1,1), grid=(GSZ,1)) 
    tac.record()
    tac.synchronize()
    time_P2P = tec.time_till(tac)*1e-3    

if mode == 'cpu':   
    tec = time.time()
    P2P_c_cy(xtHost, ytHost, ztHost, xsHost, ysHost, zsHost, intPtrHost, offsetTarHost, sizeTarHost, offsetSrcHost, offsetIntHost, paux, pxaux, pyaux, pzaux, mHost, kappa, eps)
    tac = time.time()
    time_P2P = tac - tec


if mode == 'pycuda':
    tec.record()
    pDev.get(paux)
    pxDev.get(pxaux)
    pyDev.get(pyaux)
    pzDev.get(pzaux)
    tac.record()
    tac.synchronize()
    time_trans += tec.time_till(tac)*1e-3    

if mode == 'pycuda':
    tec.record()
if mode == 'cpu':
    tec = time.time()

t = -1
for CI in twig:
    t += 1
    CI_start = offsetTarHost[t]
    CI_end = offsetTarHost[t] + sizeTarHost[t]
    targets = tarPtr[CI_start:CI_end]
    p[targets] += paux[CI_start:CI_end]
    px[targets] += pxaux[CI_start:CI_end]
    py[targets] += pyaux[CI_start:CI_end]
    pz[targets] += pzaux[CI_start:CI_end]

if mode == 'pycuda':
    tac.record()
    tac.synchronize()
    time_unsort = tec.time_till(tac)*1e-3    
if mode == 'cpu':
    tac = time.time()
    time_unsort = tac - tec

toc = time.time()
time_eval = toc-tic
    
end_tree = time.time()
print 'Mode: '+mode
print 'error phi  : %s'%sqrt(sum((phi_d[0:sample_size]-p[0:sample_size])**2)/sum(phi_d[0:sample_size]**2))
print 'error phi_x: %s'%sqrt(sum((phi_xd[0:sample_size]-px[0:sample_size])**2)/sum(phi_xd[0:sample_size]**2))
print 'error phi_y: %s'%sqrt(sum((phi_yd[0:sample_size]-py[0:sample_size])**2)/sum(phi_yd[0:sample_size]**2))
print 'error phi_z: %s\n'%sqrt(sum((phi_zd[0:sample_size]-pz[0:sample_size])**2)/sum(phi_zd[0:sample_size]**2))
print 'Time Direct  : %f\n'%time_direct
print 'Time tree gen: %f'%time_tree
print 'Time indices : %f'%time_ind
print 'Time sources : %f'%time_src
print 'Time precomp : %f'%time_pre
print 'Time P2M     : %f'%time_P2M
print 'Time M2M     : %f'%time_M2M
print 'Time packing : %f'%time_pack
print 'Time evaluate: %f'%time_eval
print '   Data trans: %f'%time_trans
print '   Time M2P  : %f'%time_M2P
print '   Time P2P  : %f'%time_P2P
print '-----------------------'
print 'TOTAL TIME   : %f'%(end_tree-start_tree)
