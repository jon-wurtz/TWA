# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 11:13:35 2016

@author: Jonathan Wurtz
"""
from hamiltonian import *
from fabc import *
from numpy import *
import sys

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
        
    def product_IC(self,psi):
        '''
        Creates a product state IC, where psi is either a vector eg for SU(3):
           a|+> + b|0> + c|-> goes to [['z',0,a],['z',1,b],['z',2,c]]
           where 'z' and the second number denote the basis vector.           
        '''
        
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
        '''
        self.mask = mask
        self.saveit = num_samples>0
        if self.saveit>0:
            self.data = zeros(num_samples)
            self.T = zeros(num_samples)
            self.index = 0
            if mask=='casmir':
                self.data = zeros([50,50,num_samples])
    
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
                
            if self.saveit:
                self.data[self.index]=X1**2 + X2**2 - XX
                self.T[self.index]=T
                self.index+=1
            return X1**2 + X2**2 - XX
        
        if self.mask=='Sz':
            if self.saveit:
                self.data[self.index]=average(data[:,:,2])
                self.T[self.index]=T
                self.index+=1
            return average(data[:,:,2])
        if self.mask=='casmir':
            if self.saveit:
                self.data[0:data.shape[0],0:data.shape[1],self.index] = sum(data**2,2)
                self.T[self.index]=T
                self.index+=1
            return sum(sum(sum(data**2,2)))
    
    def plot(self):
        if self.mask=='casmir':
            plot(self.T[0:self.index],sum(sum(sum(self.data[:,:,0:self.index]))))
        plot(self.T[0:self.index],real(self.data[0:self.index]))



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
    a nonlinear model H = n(n-1) + (a_+) + (a_-) in SU(4) representation
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


# Dimension, Size, n_iterations, IC, T,nsteps,nobs
if __name__=="__main__":
    if len(sys.argv)==1: # Run one batch with whatever settings
        params = Hubbard_SU3(2,2,-0.5,1)
        params['verbose']='f'
        di =doIT(params)
        
        T = 20.
        tstep = 0.001
        print 'Nsteps:', T/tstep
        n_samples = 1000
        bulkdata = zeros(2000)
        for i in range(10):
            print i
            di.obs.reset()
            di.product_IC([['z',1,1]])
            di.run(T,tstep,T/n_samples)
            bulkdata += di.obs.data
        di.obs.plot()
        #di.obs.data[0:di.obs.index].tofile('out.dat')
        #bulkdat[0:di.obs.index] += di.obs.data[0:di.obs.index]
    else:
        #bulkdat = zeros(2000)
        params0 = Hubbard_SU3(int(sys.argv[1]),int(sys.argv[2]),-1,1)
            
        print 'Bose-Hubbard Model: J, U, dim, n terms,IC,n points, Tstep'
        print -params0['terms'][0][0][1]
        print params0['terms'][-1][0][1]*(-6/sqrt(3))
        print params0['dim'][0]
        print params0['dim'][1]
        print sys.argv[4]
        print int(sys.argv[3])
        print double(sys.argv[5])/double(sys.argv[6])
        
        for i in range(int(sys.argv[3])):
            di = doIT(params0)
            if sys.argv[4]=='mott':
                di.product_IC([['x',0,1]])
            elif sys.argv[4]=='super':
                di.product_IC([['z',0,1]])
            dt = double(sys.argv[5])/double(sys.argv[6])
            dtobs = double(sys.argv[5])/double(sys.argv[7])
            di.run(double(sys.argv[5]),dt,dtobs)
            #bulkdat[0:di.obs.index] += di.obs.data[0:di.obs.index]
        '''
        if len(sys.argv)==9:
            f = open(sys.argv[8],'w')
            for i in range(di.obs.index-1):
                print >>f,di.obs.T[i],str('\t'),bulkdat[i]/(i+1)
            f.close()
        '''
        
