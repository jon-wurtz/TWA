# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 11:13:35 2016

@author: Jonathan Wurtz
"""
from hamiltonian import *
from fabc import *
from numpy import *
import sys
#from numba import jit

class doIT():
    def __init__(self,params):
        '''
        params is a dict with elements:
        SU:    SU(N), the group that we're closed under
        terms: a list of terms in the Hamiltonian
        dim:   an array [dimension, num elements per side]
        obs:   from class observable(), which returns a single number,
                eg, magnitization etc. from data.
        ics:   initial conditions for the run. These can also be set
                from doIT.product_IC()                
        
        '''
        self.H = hamiltonian(params['SU'],params['dim'][0])
        #self.f_abc = fabc.SU_sparse(params['SU'])
        if params['SU']==3:
            self.f_abc = old_SU3_fabc()
            self.basis = old_SU3()
        elif params['SU']==4:
            self.f_abc = SU_4_basis()[1]
            self.basis = SU_4_basis()[0]
        else:
            raise 'SU != {3,4}'
        
        self.obs = params['obs']
        self.dim = params['dim'][0]
        
        for term in params['terms']:
            self.H.add_term(term)
        
        dims = list(ones(params['dim'][0])*params['dim'][1])
        dims.append(params['SU']**2-1)
        if not 'ics' in params:
            self.data = zeros(dims)
        else:
            if params['ics'].shape==tuple(dims):
                self.data = params['ics']
            else:
                raise 'ICs are the wrong shape!!'
            
        self.params = params
        
    def product_IC(self,psi,novariance=False):
        '''
        Creates a product state IC, where psi is either a vector eg for SU(3):
           a|+> + b|0> + c|-> goes to [['z',0,a],['z',1,b],['z',2,c]]
           where 'z' and the second number denote the basis vector.           
        '''
        
        if type(psi)==ndarray:
            if (psi.shape!=self.data.shape[0:-1]):
                raise 'Incorrect Shape!'
                
            sies = psi.shape[0]
            if self.dim==2:
                itrs = list(array(meshgrid(range(sies),range(sies))).reshape(2,sies**2).transpose())
            elif self.dim==3:
                itrs = list(array(meshgrid(range(sies),range(sies),range(sies))).reshape(3,sies**3).transpose())
            for element in itrs:
                wm = {'x':0,'y':1,'z':2}
                psi2 = zeros(self.params['SU'])*1j
                for x in psi[tuple(element)]:
                    psi2 += x[2]*linalg.eig(self.basis[:,:,wm[x[0]]])[1][:,x[1]]
                
                # Renormalize to 1...
                psi2 = psi2*(dot(psi2,conj(psi2))**-.5)
                #print psi2
                
                # Find the ICs in the SU(N) basis
                T_init = find_ICs(psi2)
                if novariance:
                    self.data[tuple(element)] = T_init[0]#random.multivariate_normal(T_init[0],0*T_init[1])
                else:
                    self.data[tuple(element)] = random.multivariate_normal(T_init[0],T_init[1])
        
        else:
            # Convert ICs to the Z-basis
            wm = {'x':0,'y':1,'z':2}
            psi2 = zeros(self.params['SU'])*1j
            for x in psi:
                psi2 += x[2]*linalg.eig(self.basis[:,:,wm[x[0]]])[1][:,x[1]]
            
            # Renormalize to 1...
            psi2 = psi2*(dot(psi2,conj(psi2))**-.5)
            #print psi2
            
            # Find the ICs in the SU(N) basis
            T_init = find_ICs(psi2)
            #print 'Determinant of Covariance:',linalg.det(T_init[1])
            if novariance:
                self.data = random.multivariate_normal(T_init[0],0*T_init[1],size=self.data.shape[0:-1])
            else:
                self.data = random.multivariate_normal(T_init[0],T_init[1],size=self.data.shape[0:-1])
    


    def run(self,T,dt,dt_obs=Inf):
        if dt>(dt_obs/2):
            print 'Its going to be slow going with this many observations'
        tt = 0
        for kk in range(int(T/dt)):
            
            if abs(dt*kk-tt)<dt_obs:
                printit = self.obs.get(self.data,tt+dt_obs)
                tt+=dt_obs
                if 't' in self.params['verbose']:
                    print tt,str('\t'),T,str('\t'),printit
            
            if self.dim==2:
              for abc in self.f_abc: # Use sparse f_abc...
                self.data[:,:,abc[0]] += dt*abc[3]*self.data[:,:,abc[2]]*self.H.dH(self.data,abc[1])
            elif self.dim==3:
              for abc in self.f_abc: # Use sparse f_abc...
                self.data[:,:,:,abc[0]] += dt*abc[3]*self.data[:,:,:,abc[2]]*self.H.dH(self.data,abc[1])
            elif self.dim==1:
              for abc in self.f_abc: # Use sparse f_abc...
                self.data[:,abc[0]] += dt*abc[3]*self.data[:,abc[2]]*self.H.dH(self.data,abc[1])
            else:
                raise 'Something has gone wrong!'
            


class observable():
    def __init__(self,mask,num_samples):
        '''
        mask: a string describing the mask
         super: Superfluid density, from Bose-Hubbard model
         casmir: the Casmir operator sum(X**2)
         Sz: the Sz observable...
         num_particles: The number of particles per state
        '''
        self.mask = mask
        self.saveit = num_samples>0
        if self.saveit>0:
            self.data = zeros(num_samples)
            self.T = zeros(num_samples)
            self.index = 0
            if mask=='casmir':
                self.data = zeros([50,50,num_samples])
            if mask=='num_particles':
                self.data = zeros([num_samples,4])
                self.SU4 = [to_SUbasis(diag([1,0,0,0])),to_SUbasis(diag([0,1,0,0])),to_SUbasis(diag([0,0,1,0])),to_SUbasis(diag([0,0,0,1]))]
                self.SU3 = [to_SUbasis(diag([1,0,0])),to_SUbasis(diag([0,1,0])),to_SUbasis(diag([0,0,1]))]
    
    def reset(self):
        if self.saveit>0:
            self.data = self.data*0
            self.T = self.data*0
            self.index = 0
            if self.mask=='casmir':
                self.data = zeros([50,50,num_samples])
    def get(self,data,T):
        if self.mask=='superfluid':
            
            if len(data.shape)==3:
                X1 = average(data[:,:,0])
                X2 = average(data[:,:,1])
                XX = average(data[:,:,0]**2 + data[:,:,1]**2)/(data.shape[0]*data.shape[1])
            elif len(data.shape)==4:
                X1 = average(data[:,:,:,0])
                X2 = average(data[:,:,:,1])
                XX = average(data[:,:,:,0]**2 + data[:,:,:,1]**2)/(data.shape[0]*data.shape[1]*data.shape[2])
            else:
                raise 'Something went Wrong!'
            
            # SU(4) vs SU(3) has a different 
            mod = 0           
            if data.shape[-1]==8:
                mod=1
            elif data.shape[-1]==15:
                mod=2.5
            else:
                raise 'BAD!'
            if self.saveit:
                self.data[self.index]=(X1**2 + X2**2 - XX)*mod
                self.T[self.index]=T
                self.index+=1
            return (X1**2 + X2**2 - XX)*mod
        
        elif self.mask=='Sz':
            if self.saveit:
                self.data[self.index]=average(data[:,:,2])
                self.T[self.index]=T
                self.index+=1
            return average(data[:,:,2])
        elif self.mask=='casmir':
            if self.saveit:
                self.data[0:data.shape[0],0:data.shape[1],self.index] = sum(data**2,2)
                self.T[self.index]=T
                self.index+=1
            return sum(sum(sum(data**2,2)))
        
        elif self.mask=='num_particles':
            if self.saveit:
                for i in range(int(sqrt(data.shape[2]+1))):
                    mm = zeros(int(sqrt(data.shape[2]+1)))
                    mm[i] = 1
                    self.data[self.index,i] = average(average(dot(data,to_SUbasis(diag(mm)))))
                    self.T[self.index]=T
                self.index+=1
            return None
        
        elif self.mask=='Sz2':
            if self.saveit:
                if data.shape[2]==8:
                    self.data[self.index]=average(data[:,:,7])*(-0.57735027)+2/3.
                    self.T[self.index]=T
                    self.index+=1
                    return average(data[:,:,7])*(-0.57735027)+2/3.
                elif data.shape[2]==15:
                    self.data[self.index]=average(data[:,:,7])*sqrt(2)
                    self.T[self.index]=T
                    self.index+=1
                    return average(data[:,:,7])*sqrt(2)
        
        elif self.mask=='diffusion':
            if data.shape[2]==8:
                self.data[self.index] = (data[0,0,2]+1)/sum(data[:,:,2]+1)
            if data.shape[2]==15:
                self.data[self.index] = (data[0,0,2]*1.5811383+1.5)/sum(data[:,:,2]*1.5811383+1.5)
                
            self.T[self.index]=T
            self.index+=1
            return self.data[self.index-1]
        else:
            raise 'Bad mask'

                
    
    def plot(self):
        if self.mask=='casmir':
            plot(self.T[0:self.index],sum(sum(sum(self.data[:,:,0:self.index]))))
        plot(self.T[0:self.index],real(self.data[0:self.index]))

    def put(self,filename):
        f = open(filename,'a')
        print >>f,'---DATA---'
        print >>f,self.data.tostring()
        f.close()


def Hubbard_SU3(dim,sies,J,U):
    output = {}
    output['SU']=3
    output['dim']=[dim,sies]
    output['obs']=observable('superfluid',2000)
    output['verbose']='t'
    
    terms = []
    
    
    if output['dim'][0]==3:
        terms.append([[0,J],[0,array([0,0,1])]])
        terms.append([[0,J],[0,array([0,1,0])]])
        terms.append([[0,J],[0,array([1,0,0])]])    
        
        terms.append([[1,J],[1,array([0,0,1])]])
        terms.append([[1,J],[1,array([0,1,0])]])
        terms.append([[1,J],[1,array([1,0,0])]])
        
    elif output['dim'][0]==2:
        terms.append([[0,J],[0,array([0,1])]])
        terms.append([[0,J],[0,array([1,0])]])
        
        terms.append([[1,J],[1,array([0,1])]])
        terms.append([[1,J],[1,array([1,0])]])
        
    
    terms.append([[7,-U*sqrt(3)/6]])
    
    
    output['terms']=terms
    return output

def local_Hubbard_SU4():
    '''
    a nonlinear model H = n(n-2) + (a_+) + (a_-) in SU(4) representation
    '''
    output = {}
    output['SU'] = 4
    output['dim'] = [2,40] # Although local, having finite size
                           # lets us gather statistics all in 1 run.
    output['obs'] = observable('Sz',10000)
    output['verbose']='t'
    
    terms = []
    
    # Sz^2 + Sz
    terms.append([[2,sqrt(5/2.)]])
    terms.append([[7,sqrt(2.)]])
    
    # Sx
    terms.append([[0,sqrt(5/2.)]])

    output['terms'] = terms
    return output

def Hubbard_SUN(dim,sies,J,U,N,meanfield=False,mu_fname=None):
    '''
    Paramaters for a SU(N) Bose-Hubbard model.
    here, we use the matrices for boson raising and lowering operators,
     and NOT those for spin raising and lowering.
    
    dim - dimension: in {1,2,3}
    sies - size of lattice. Please, <30?
    J - Coupling strength. Should be (-)
    U - Interaction strenth. Should be (+)
    N - SU(*). in {3,4}
        - Can also be:
          - "old" - SU(3) rep of Spin operators, not bosons
          - "SU3" - SU3 Hamiltonian embedded in SU4 (w/ bosons)
    
    Note that the zeroth element of the wavefunc'n is
     the fully occupied mode, and the last element is the unoccupied one.
    
    '''
    
    
    output = {}
    if N=='old':
        output['SU'] = 3
    elif N=='SU3':
        output['SU'] = 4
    else:
        output['SU'] = N
    output['dim'] = [dim,sies]
    output['obs'] = observable('superfluid',10000)
    output['verbose']='f'
    
    terms = []
    
    # Define matrices...
    if N==4:
        matr_a = diag([sqrt(3),sqrt(2),1],k=1) # the a+ operator
        matr_U = diag([3,0,-1,0]) # U*n(n-2)
        matr_mu = diag([3,2,1,0])
    elif N==3:
        matr_a = diag([sqrt(2),1],k=1) # the a+ operator
        matr_U = diag([0,-1,0]) # U*n(n-1)
        matr_mu = diag([2,1,0])
    elif N=='old':
        matr_a = diag([sqrt(2),sqrt(2)],k=1) # the S+ operator
        matr_U = diag([0,-1,0]) # U*n(n-2)
    elif N=='SU3': # Use SU4 generators but SU3 Hamiltonian
        matr_a = diag([0,sqrt(2),1],k=1) # the a+ operator but with no 4th element.
        matr_U = diag([0,0,-1,0]) # U*n(n-2)
    else:        
        raise 'Bad SU(N)! N={3,4,old,SU3}'
    
    # Define nonlocal interactions
    UU = to_SUbasis(matr_a)
    VV = to_SUbasis(matr_a.transpose())
    UV = real(outer(UU,VV) + outer(VV,UU)) # Assert: it is real.
    
    
    neighbors = list(identity(dim).astype(int32)) # Define nearest-neighbor interaction
   
    if meanfield:
        print 'Doing it with Mean Field!!'
        # Mean-Field!
        J =  0.25*J/(sies**dim-1) # Everything is quadruple-counted; once by inversion symmetry and once by i<->j
        if dim==2:
            neighbors = list(array((meshgrid(range(sies),range(sies)))).reshape(2,sies**2).transpose())[1:-1]
        elif dim==3:
            neighbors = list(array((meshgrid(range(sies),range(sies),range(sies)))).reshape(2,sies**3).transpose())[1:-1]
    
    if J!=0:
        for neighbor in neighbors:
            for UV_ in list(array(nonzero(UV)).transpose()): # I'm an abomination
                terms.append([[int(UV_[0]),J*UV[UV_[0],UV_[1]]],[int(UV_[1]),neighbor]])
        
    # Define local interactions...
    Q = real(to_SUbasis(matr_U)) # Assert: it is real.
    for Q_ in list(array(nonzero(Q)).transpose()):
        terms.append([ [int(Q_[0]), U*0.5*Q[Q_[0]]] ])
    
    # Define Chemical Potential (if its there)
    if mu_fname!=None and type(N)==int:
        f = open(mu_fname,'rb')
        print 'Loading chemical potential data!'
        mu = loadtxt(f)
        f.close()
        R = real(to_SUbasis(matr_mu))
        for R_ in list(array(nonzero(R)).transpose()):
            terms.append([ [int(R_[0]), mu*R[R_[0]]] ])
    
    output['terms'] = terms
    return output

        
