# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 11:08:02 2016

@author: Jonathan Wurtz
"""
from numpy import *
# Defines the structure constants for SU(N) in the fundemental representation

def old_SU3():
    '''
    Basis matrices for SU(3)
    '''
    M = zeros([3,3,8])*1j
    M[:,:,0] = 2**-.5*array([[0,1,0],[1,0,1],[0,1,0]])
    M[:,:,1] = 1j*2**-.5*array([[0,-1,0],[1,0,-1],[0,1,0]])
    M[:,:,2] = array([[1,0,0],[0,0,0],[0,0,-1]])
    M[:,:,3] = array([[0,0,1],[0,0,0],[1,0,0]])
    M[:,:,4] = 1j*array([[0,0,-1],[0,0,0],[1,0,0]])
    M[:,:,5] = 2**-.5*array([[0,-1,0],[-1,0,1],[0,1,0]])
    M[:,:,6] = 1j*2**-.5*array([[0,1,0],[-1,0,-1],[0,1,0]])
    M[:,:,7] = 3**-.5*array([[-1,0,0],[0,2,0],[0,0,-1]])
    #M[:,:,8] = array([[1,0,0],[0,1,0],[0,0,1]]) # Identity...
    
    return M

def old_SU3_fabc():
    '''
    Structure functions for SU(3)
    '''
    f_abc = []
    # Non-zero structure constants...
    f_gen = []
    f_gen.append([1,2,3,1])
    f_gen.append([1,4,7,1])
    f_gen.append([1,6,5,1])
    f_gen.append([2,4,6,1])
    f_gen.append([2,5,7,1])
    f_gen.append([3,6,7,1])
    f_gen.append([1,7,8,sqrt(3)])
    f_gen.append([2,8,6,sqrt(3)])
    f_gen.append([3,4,5,2])
    
    # Permute 'em here...
    for gen in f_gen:
        for i in range(3):
            f_abc.append([gen[(0+i)%3]-1,gen[(1+i)%3]-1,gen[(2+i)%3]-1,gen[3]])
            f_abc.append([gen[(0+i)%3]-1,gen[(2+i)%3]-1,gen[(1+i)%3]-1,-gen[3]])
    
    return f_abc

def find_ICs(ICs):
    '''
    ICs is a vector of length 3 or 4, are the normalized
    variables in the wavefunction.
    It returns a vector and matrix, where the vector is the center
    of a gaussian distribution in the Weyl-transformed SU(N) generator
    phase space coordinates, and the matrix is the covarience of
    a gaussian distribution, or equivilently the correlations and
    fluctuations between variables.
    '''
    if len(ICs)==3:
        out = zeros(8)*1j
        M = old_SU3()
    elif len(ICs)==4:
        out = zeros(15)*1j
        M = SU_4_basis()[0]
    else:
        raise 'SU !={3,4}'
        
    for i in range(len(out)):
        out[i] = dot(conj(ICs),dot(M[:,:,i],ICs))
        
    var = zeros([len(out),len(out)])*1j
    for i in range(len(out)):
        for j in range(len(out)):
            var[i,j] = 0.5*dot(conj(ICs),dot(dot(M[:,:,i],M[:,:,j])+dot(M[:,:,j],M[:,:,i]),ICs))-out[i]*out[j]
            if abs(imag(var[i,j]))>1e-10:
                raise 'Bad variance!'
    
    out = real(out)
    var = real(var)
    out[abs(out)<1e-10]=0
    var[abs(var)<1e-10]=0
    return out,var

def SU_4_basis(checkit=False):
    '''
    Based off of Kiselev Mikhail's work
    "Controllable 1=N expansion for SU(N) quantum dynamical problems"
    
    Returns:
    - the basis matrices as described by KM as a structure M[:,:,i] where
       i is the ith basis matrix
    - the sparse structure constants associated with the particular choice
       of basis
    - Note that with a quick modification (fabc -> f_abc) the structure constants
       can be made not-sparse.
      
    You can check if the trace identities, Jacoby identity, and permutations
     work by passing checkit=True.    
    '''
    N=4
    M = zeros([4,4,15])*1j
    sqt = sqrt(3)/2
    Jx = array([[0,sqt,0,0],[sqt,0,1,0],[0,1,0,sqt],[0,0,sqt,0]])
    Jy = 1j*array([[0,sqt,0,0],[-sqt,0,1,0],[0,-1,0,sqt],[0,0,-sqt,0]])
    Jz = array([[3/2.,0,0,0],[0,1/2.,0,0],[0,0,-1/2.,0],[0,0,0,-3/2.]])
    
    # Sx
    M[:,:,0] = sqrt(2./5)*Jx
    #Sy    
    M[:,:,1] = sqrt(2./5)*Jy
    #Sz    
    M[:,:,2] = sqrt(2./5)*Jz
    
    M[:,:,3] = 6**-.5*(dot(Jx,Jx)-dot(Jy,Jy))
    M[:,:,4] = 6**-.5*(dot(Jx,Jy)+dot(Jy,Jx))
    M[:,:,5] = 6**-.5*(dot(Jx,Jz)+dot(Jz,Jx))
    M[:,:,6] = 6**-.5*(dot(Jy,Jz)+dot(Jz,Jy))
    M[:,:,7] = 18**-.5*(2*dot(Jz,Jz)-dot(Jy,Jy)-dot(Jx,Jx))
    M[:,:,8] = 6**-1*(2*ac(Jx,dot(Jx,Jx)-dot(Jy,Jy)) - ac(Jz,ac(Jx,Jz)) + ac(Jx,2*dot(Jz,Jz)-dot(Jy,Jy)-dot(Jx,Jx)))
    M[:,:,9] = 6**-1*(2*ac(Jy,dot(Jx,Jx)-dot(Jy,Jy)) + ac(Jz,ac(Jy,Jz)) - ac(Jy,2*dot(Jz,Jz)-dot(Jy,Jy)-dot(Jx,Jx)))
    
    M[:,:,10] = 6**-.5*ac(Jz,dot(Jx,Jx)-dot(Jy,Jy))
    M[:,:,11] = 6**-.5*ac(Jz,ac(Jx,Jy))
    
    M[:,:,12] = 60**-.5*(ac(Jz,ac(Jx,Jz)) + ac(Jx,2*dot(Jz,Jz)-dot(Jy,Jy)-dot(Jx,Jx)))
    M[:,:,13] = 60**-.5*(ac(Jz,ac(Jy,Jz)) + ac(Jy,2*dot(Jz,Jz)-dot(Jy,Jy)-dot(Jx,Jx)))
    
    M[:,:,14] = (10./9)**.5*(dot(dot(Jz,Jz),Jz)- 41./20*Jz)

    # Find structure functions
    f_abc = zeros([N**2-1,N**2-1,N**2-1])
    fabc = []
    for k in range(N**2-1):
        for l in range(N**2-1):
            new = dot(M[:,:,k],M[:,:,l])-dot(M[:,:,l],M[:,:,k])
            new[abs(new)<1e-10]=0
            for m in range(N**2-1):
                f_abc[k,l,m] = trace(dot(M[:,:,m],new))/(2*1j)
                if abs(f_abc[k,l,m])>1e-6:
                    fabc.append([k,l,m,f_abc[k,l,m]])
                    
    
    # Below are double-checks to make sure everything is happy...
    if checkit:
        # Check trace identities...
        Tr = zeros([15,15])
        for i in range(Tr.shape[0]):
            for j in range(Tr.shape[1]):
                Tr[i,j] = trace(dot(M[:,:,i],M[:,:,j]))
        #imshow(Tr,interpolation='none')
        
        print 'Trace satisfied:',abs(sum(sum(Tr-identity(15)*2)))<1e-6
        
        # Find Structure Functions and check Jacobi Identity
        
        # Copy-Pasted from SU_N.py
        N=4
    
                
        f_abc[abs(f_abc)<1e-6]=0
        #f_abc = around(f_abc,6)
        # Recognize integers that have had funny floating point errors...
        #f_abc[around(around(f_abc,6),0)==around(f_abc,6)] = around(f_abc[around(around(f_abc,6),0)==around(f_abc,6)],6)
        
        # Check the Jacobi identity...
        # Claim: it works
        # Evidence: it works for SU(6) but is very slow because Python.
        #  - Speed scales as N^10...
        tot = zeros(ones(4)*N**2)
        for a in range(N**2-1):
            #print a,N**2
            for b in range(N**2-1):
                for c in range(N**2-1):
                    for e in range(N**2-1):
                        for d in range(N**2-1):
                            tot[a,b,c,e] += f_abc[a,b,d]*f_abc[d,c,e]
                            tot[a,b,c,e] += f_abc[b,c,d]*f_abc[d,a,e]
                            tot[a,b,c,e] += f_abc[c,a,d]*f_abc[d,b,e]
        print 'Jacoby Identity Satisfied:',(abs(tot)<1e-3).all()
        
        
        # Check if things are symmetric under permutations...
        truth = True
        for i in range(f_abc.shape[0]):
            for j in range(f_abc.shape[1]):
                for k in range(f_abc.shape[2]):
                    truth = truth and abs(f_abc[i,j,k] - f_abc[k,i,j])<1e-7
                    truth = truth and abs(f_abc[i,j,k] - f_abc[j,k,i])<1e-7
                    truth = truth and abs(f_abc[i,j,k] - -f_abc[j,i,k])<1e-7
                    truth = truth and abs(f_abc[i,j,k] - -f_abc[k,j,i])<1e-7
        
        print 'Anti-symmetric under permutations:',truth
        
    # f_abc is also here, which is a not-sparse-matrix.
    return M,fabc
    

def to_SUbasis(matr):
    '''
    Converts a square 3x3 or 4x4 matrix into its equivilent
    vector of SU(N) operators.
    '''
    if matr.shape[0]==4:
        M = SU_4_basis()[0]
    elif matr.shape[0]==3:
        M = old_SU3()
    else:
        raise 'Bad matrix shape: N={3,4}'
    
    out = zeros(M.shape[2])*1j
    for i in range(M.shape[2]):
        out[i] = 0.5*trace(dot(M[:,:,i],matr))
    
    out[abs(out)<1e-10]=0
    return out
        
        
        
    
        
    
def ac(A,B):
    '''
    Anticommutator
    '''
    return dot(A,B)+dot(B,A)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
