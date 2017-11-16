from numpy import *
from numpy import sum as npsum
from scipy.misc import factorial
#from scipy.misc.common import comb
from scipy.misc import comb

# Wrapped code
from P2Pwrap import direct_c
from calculateMultipoles import P2M, M2M
from M2Pwrap import getCoeff, setIndex, setIndex_arr

import time

class Cell():
    def __init__ (self, NCRIT, Nm):
        self.nsource = 0       # Number of source particles
        self.ntarget = 0       # Number of target particles
        self.nchild  = 0      # Number of child boxes in binary
                              # This will be a 8bit value and if certain 
                              # child exists, that bit will be 1.

        self.source = []                            # Pointer to source particles
        self.target = []                            # Pointer to target particles
        self.xc = 0.                                # x position of cell
        self.yc = 0.                                # y position of cell
        self.zc = 0.                                # z position of cell
        self.r  = 0.                                # cell radius

        self.parent = 0                             # Pointer to parent cell
        self.child  = zeros(8, dtype=int32)         # Pointer to child cell

        self.M = zeros(Nm)                          # Array with multipoles

        self.P2P_list = array([], dtype=int32)       # Pointer to cells that interact with P2P
        self.M2P_list = array([], dtype=int32)       # Pointer to cells that interact with M2P
        self.list_ready = 0                          # Flag to know if P2P list is already generated

def add_child(octant, Cells, i, NCRIT, Nm):
    # add_child adds child cell to Cells array
    # octant: octant of the child cell
    # Cells : arrays with cells
    # i     : index of parent cell in Cells array

    CN    = Cell(NCRIT, Nm) # CN: child cell
    CN.r  = Cells[i].r/2
    CN.xc = Cells[i].xc + CN.r*((octant&1)*2-1) # octant&X returns X if true
    CN.yc = Cells[i].yc + CN.r*((octant&2)-1)   # Want to make ((octant&X)*Y - Z)=1
    CN.zc = Cells[i].zc + CN.r*((octant&4)/2-1)
    CN.parent = i
    Cells[i].child[octant] = len(Cells)
    Cells[i].nchild|=(1<<octant)
    Cells.append(CN)

def split_cell(x, y, z, Cells, C, NCRIT, Nm):
    # split_cell splits cell with more than NCRIT particles
    # x,y,z: positions of particles
    # Cells: array of cells
    # C    : index of cell to be split in Cells array

    for l in Cells[C].target:
        octant = (x[l]>Cells[C].xc) + ((y[l]>Cells[C].yc) << 1) + ((z[l]>Cells[C].zc) << 2)
        if (not(Cells[C].nchild & (1<<octant))): # Ask if octant exists already
            add_child(octant, Cells, C, NCRIT, Nm)

        CC = Cells[C].child[octant] # Pointer to child cell
        Cells[CC].target.append(l)
        Cells[CC].ntarget += 1

        if (Cells[CC].ntarget >= NCRIT):
            split_cell(x, y, z, Cells, CC, NCRIT, Nm)

def generateTree(xi, yi, zi, NCRIT, Nm, N, x0):

    C0 = Cell(NCRIT, Nm)
    C0.xc = x0[0]
    C0.yc = x0[1]
    C0.zc = x0[2]
    C0.r  = x0[3]

    Cells = []
    Cells.append(C0)

    for i in range(N):

        C = 0 
        while (Cells[C].ntarget>=NCRIT):
            Cells[C].ntarget+=1

            octant = int(xi[i]>Cells[C].xc) + int(yi[i]>Cells[C].yc)*2 + int(zi[i]>Cells[C].zc)*4
            if (not(Cells[C].nchild & (1<<octant))):
                add_child(octant, Cells, C, NCRIT, Nm)
        
            C = Cells[C].child[octant]

        Cells[C].target.append(i) 
        Cells[C].ntarget += 1

        if (Cells[C].ntarget>=NCRIT):
            split_cell(xi,yi,zi,Cells,C, NCRIT, Nm)

    return Cells

def computeIndices(P):
    II = []
    JJ = []
    KK = []
    index = []
    index_large = zeros((P+1)*(P+1)*(P+1), dtype=int32)
    for ii in range(P+1):
        for jj in range(P+1-ii):
            for kk in range(P+1-ii-jj):
                index.append(setIndex(P,ii,jj,kk))
                II.append(ii)
                JJ.append(jj)
                KK.append(kk)
                index_large[(P+1)*(P+1)*ii+(P+1)*jj+kk] = index[-1]

    II = array(II,int32)
    JJ = array(JJ,int32)
    KK = array(KK,int32)
    index = array(index,int32)
#    index = setIndex_arr(P,II,JJ,KK)
    
    return II, JJ, KK, index, index_large



def getMultipole(Cells, C, x, y, z, m, twig, II, JJ, KK, index, P, NCRIT):
    # Cells     : array of cells
    # C         : index of cell in Cells array 
    # x,y,z     : position of sources
    # m         : weight of particles
    # P         : order of Taylor expansion
    # NCRIT     : max number of source particles per cell
    # II,JJ,KK  : x,y,z powers of multipole expansion
    # index     : 1D mapping of II,JJ,KK (index of multipoles)

    if (Cells[C].ntarget>=NCRIT):
        for c in range(8):
            if (Cells[C].nchild & (1<<c)):
                twig = getMultipole(Cells, Cells[C].child[c], x, y, z, m, twig, II, JJ, KK, index, P, NCRIT)
    else:
        '''
        for l in Cells[C].source:
            dx = Cells[C].xc - x[l]
            dy = Cells[C].yc - y[l]
            dz = Cells[C].zc - z[l]

            Cells[C].M[index] += m[l]*(dx**II)*(dy**JJ)*(dz**KK)
        '''
        '''
        l = Cells[C].source
        dx = transpose(ones((len(II),Cells[C].nsource))*(Cells[C].xc - x[l]))
        dy = transpose(ones((len(II),Cells[C].nsource))*(Cells[C].yc - y[l]))
        dz = transpose(ones((len(II),Cells[C].nsource))*(Cells[C].zc - z[l]))
        constant = transpose((dx**II)*(dy**JJ)*(dz**KK))

        Cells[C].M[index] += sum(m[l]*constant, axis=1)
        '''
        l = Cells[C].source
        P2M(Cells[C].M, x[l], y[l], z[l], m[l], Cells[C].xc, Cells[C].yc, Cells[C].zc, II, JJ, KK)

        twig.append(C)

    return twig

def addSources2(x,y,z,j,Cells,C,NCRIT):
    # x,y,z: location of sources
    # j    : index of source
    # Cells: array with cells
    # C    : index of cell in Cells array
    # NCRIT: max number of source particles per cell

    if (Cells[C].ntarget>=NCRIT):
        octant = int(x[j]>Cells[C].xc) + int(y[j]>Cells[C].yc)*2 + int(z[j]>Cells[C].zc)*4
        if (Cells[C].nchild & (1<<octant)): # If child cell exists, use
            O = octant
        else:                               # If child cell doesn't exist add to closest existing child
            r = []
            child = []
            for c in range(8):
                if (Cells[C].nchild & (1<<c)):
                    dx = x[j]-Cells[Cells[C].child[c]].xc
                    dy = y[j]-Cells[Cells[C].child[c]].yc
                    dz = z[j]-Cells[Cells[C].child[c]].zc
                    r.append(sqrt(dx*dx+dy*dy+dz*dz))
                    child.append(c)
            close_child = r.index(min(r))   # Find index of closest child
            O = child[close_child]              
        addSources2(x,y,z,j,Cells, Cells[C].child[O],NCRIT)
        
    else:
        Cells[C].nsource += 1
        Cells[C].source.append(j)

def addSources(x,y,z,Cells,twig):
    # x,y,z: location of targets
    # Cells: array with cells
    # twig : array with pointers to twigs of cells array

    dx = zeros((len(twig),len(x)))
    dy = zeros((len(twig),len(x)))
    dz = zeros((len(twig),len(x)))
    j = 0
    for t in twig:
        dx[j] = x - Cells[t].xc
        dy[j] = y - Cells[t].yc
        dz[j] = z - Cells[t].zc
        j+=1
    r = sqrt(dx*dx+dy*dy+dz*dz)

    close_twig = argmin(r,axis=0)

    for j in range(len(close_twig)):
        Cells[twig[close_twig[j]]].nsource += 1
        Cells[twig[close_twig[j]]].source.append(j)


def precompute_terms(P, II, JJ, KK):
    # Precompute terms for 
    combII = array([],dtype=int32)
    combJJ = array([],dtype=int32)
    combKK = array([],dtype=int32)
    IImii  = array([],dtype=int32)
    JJmjj  = array([],dtype=int32)
    KKmkk  = array([],dtype=int32)
    index  = array([], dtype=int32)
    index_ptr = zeros(len(II)+1, dtype=int32)
    for i in range(len(II)):
        ii,jj,kk = mgrid[0:II[i]+1:1,0:JJ[i]+1:1,0:KK[i]+1:1].astype(int32)
        ii,jj,kk = ii.ravel(), jj.ravel(), kk.ravel()
        index_aux = zeros(len(ii), int32)
        setIndex_arr(P, len(ii), index_aux, ii, jj, kk)
        index = append(index, index_aux)
        index_ptr[i+1] = len(index_aux)+index_ptr[i]
        combII = append(combII, comb(II[i],ii))
        combJJ = append(combJJ, comb(JJ[i],jj))
        combKK = append(combKK, comb(KK[i],kk))
        IImii = append(IImii, II[i]-ii)
        JJmjj = append(JJmjj, JJ[i]-jj)
        KKmkk = append(KKmkk, KK[i]-kk)
    return combII, combJJ, combKK, IImii, JJmjj, KKmkk, index, index_ptr

def upwardSweep(Cells, CC, PC, P, II, JJ, KK, Index, combII, combJJ, combKK, IImii, JJmjj, KKmkk, index, index_ptr):
    # Cells     : array of cells
    # CC        : index of child cell in Cells array
    # PC        : index of parent cell in Cells array
    # P         : order of Taylor expansion
    # II,JJ,KK  : x,y,z powers of multipole expansion
    # index     : 1D mapping of II,JJ,KK (index of multipoles)

    dx = Cells[PC].xc - Cells[CC].xc
    dy = Cells[PC].yc - Cells[CC].yc
    dz = Cells[PC].zc - Cells[CC].zc

    M2M(Cells[PC].M, Cells[CC].M, dx, dy, dz, II, JJ, KK, combII, combJJ, combKK, IImii, JJmjj, KKmkk, index, index_ptr)
    ''' 
    for ij in range(len(Index)):
        ptr = arange(index_ptr[ij],index_ptr[ij+1])
        Cells[PC].M[ij] += npsum(Cells[CC].M[index[ptr]]*combII[ptr]*combJJ[ptr]*combKK[ptr]*dx**IImii[ptr]*dy**JJmjj[ptr]*dz**KKmkk[ptr])
    '''
   
def packData(Cells, CJ, CI, intPtr, mltPtr, offInt, offMlt, theta, NCRIT):
    # Cells     : array of Cells
    # CJ        : index of source cell
    # CI        : index of target cell
    # srcPtr    : array with pointers to sources 
    # mltPtr    : array with pointers to cells for M2P 
    # offSrc    : array with offsets to sources
    # offMlt    : array with offsets to cells for M2P
    # theta     : MAC criteron 
    # NCRIT     : max number of particles per cell

    if (Cells[CJ].ntarget>=NCRIT):
        for c in range(8):
            if (Cells[CJ].nchild & (1<<c)):
                CC = Cells[CJ].child[c]  # Points at child cell
                dxi = Cells[CC].xc - Cells[CI].xc
                dyi = Cells[CC].yc - Cells[CI].yc
                dzi = Cells[CC].zc - Cells[CI].zc
                r   = sqrt(dxi*dxi+dyi*dyi+dzi*dzi)
                if Cells[CI].r+Cells[CC].r > theta*r: # Max distance between particles
                    intPtr, mltPtr, offInt, offMlt = packData(Cells, CC, CI, intPtr, mltPtr, offInt, offMlt, theta, NCRIT)
                else:
                    mltPtr[offMlt] = CC
                    offMlt += 1

    else: # Else on a twig cell
        intPtr[offInt] = CJ
        offInt += 1

    return intPtr, mltPtr, offInt, offMlt
