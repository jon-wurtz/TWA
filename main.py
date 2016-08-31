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
                from doIT.wigner_IC()                
        
        '''
        self.H = hamiltonian(params['SU'],params['dim'][0])
        #self.f_abc = fabc.SU_sparse(params['SU'])
        self.f_abc = old_SU3_fabc()
        self.obs = params['obs']
        self.dim = params['dim'][0]
        
        for term in params['terms']:
            self.H.add_term(term)
        self.H.print_hamiltonian()
        
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
    def wigner_IC(self,ICs):
        if ICs.shape!=self.data.shape:
            raise 'Wrong ICs shape!'
        self.data = ICs
        

    def run(self,T,dt,dt_obs=Inf):
        if dt>(dt_obs/2):
            print 'Its going to be slow going with this many observations'
        tt = 0
        for kk in range(int(T/dt)):
            
            if abs(dt*kk-tt)<dt_obs:
                self.obs.get(self.data,tt+dt_obs)
                tt+=dt_obs
                if 't' in self.params['verbose']:
                    print tt,str('\t'),T,str('\t'),self.obs.data[self.obs.index-1]
            
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
        '''
        self.mask = mask
        self.data = zeros(num_samples)
        self.T = zeros(num_samples)
        self.index = 0
        if mask=='casmir':
            self.data = zeros([50,50,num_samples])
    def get(self,data,T):
        if self.mask=='superfluid':
            if len(data.shape)==3:
                X1 = average(data[:,:,0])
                X2 = average(data[:,:,1])
            elif len(data.shape)==4:
                X1 = average(data[:,:,:,0])
                X2 = average(data[:,:,:,1])
            else:
                raise 'Something went Wrong!'
            self.data[self.index]=X1**2 + X2**2
            self.T[self.index]=T
            self.index+=1
        
        if self.mask=='Sz2':
            self.data[self.index]=average(data[:,:,0])
            self.T[self.index]=T
            self.index+=1
        if self.mask=='casmir':
            self.data[0:data.shape[0],0:data.shape[1],self.index] = sum(data**2,2)
            self.T[self.index]=T
            self.index+=1
    
    def plot(self):
        plot(self.T[0:self.index],real(self.data[0:self.index]))



def Hubbard_old(dim,sies,inits):
    output = {}
    output['SU']=3
    output['dim']=[dim,sies]
    output['obs']=observable('superfluid',2000)
    terms = []
    if inits=='super':
        J= -1.0
        U = 1
    elif inits=='mott':
        J = -.0025
        U= 1
    else:
        raise 'Bad IC'
    
    #IC[:,:,:,7] = 1/sqrt(3) 
    
    
    if output['dim'][0]==3:
        terms.append([[0,J],[0,array([0,0,1])]])
        terms.append([[0,J],[0,array([0,1,0])]])
        terms.append([[0,J],[0,array([1,0,0])]])    
        
        terms.append([[1,J],[1,array([0,0,1])]])
        terms.append([[1,J],[1,array([0,1,0])]])
        terms.append([[1,J],[1,array([1,0,0])]])
        IC = zeros([output['dim'][1],output['dim'][1],output['dim'][1],8])

    elif output['dim'][0]==2:
        terms.append([[0,J],[0,array([0,1])]])
        terms.append([[0,J],[0,array([1,0])]])
        
        terms.append([[1,J],[1,array([0,1])]])
        terms.append([[1,J],[1,array([1,0])]])
        
        IC = zeros([output['dim'][1],output['dim'][1],8])
        
    
    terms.append([[7,-U*sqrt(3)/6]])
    
    
    M = old_SU3()
    if inits=='mott':
        T_init = find_ICs(linalg.eig(M[:,:,2])[1][:,1]) # |Z=0>
        IC += random.multivariate_normal(T_init[0],T_init[1],size=IC.shape[0:-1])
    elif inits=='super':
        T_init = find_ICs(linalg.eig(M[:,:,0])[1][:,0]) # |X=1>
        IC += random.multivariate_normal(T_init[0],T_init[1],size=IC.shape[0:-1])
        
    output['ics'] = IC
    
    output['verbose']='t'    
    
    output['terms']=terms
    return output

    
    


# Dimension, Size, n_iterations, IC, T,nsteps,nobs, fname
if __name__=="__main__":
    if len(sys.argv)==1: # Run one batch with whatever settings
        di =doIT(Hubbard_old(3,10,'super'))
        T = 40.
        nsteps = 10000
        n_samples = 1000
        di.run(T,T/nsteps,T/n_samples)
        di.obs.plot()
        #di.obs.data[0:di.obs.index].tofile('out.dat')
        #bulkdat[0:di.obs.index] += di.obs.data[0:di.obs.index]
    else:
        #bulkdat = zeros(2000)
        params0 = Hubbard_old(int(sys.argv[1]),int(sys.argv[2]),sys.argv[4])
        print 'Bose-Hubbard Model: J, U, dim, n terms,IC,n points, Tstep'
        print -params0['terms'][0][0][1]
        print params0['terms'][-1][0][1]*(-6/sqrt(3))
        print params0['dim'][0]
        print params0['dim'][1]
        print sys.argv[4]
        print int(sys.argv[3])
        print double(sys.argv[5])/double(sys.argv[6])
        
        for i in range(int(sys.argv[3])):
            di = doIT(Hubbard_old(int(sys.argv[1]),int(sys.argv[2]),sys.argv[4]))
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
        
