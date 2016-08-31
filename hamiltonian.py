# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 10:49:46 2016

@author: Jonathan Wurtz
"""

# Here, we build the class that defines the total Hamiltonian
#  of the system.
# We assume translational symmetry for lattice interactions.
from numpy import *

class hamiltonian():
    def __init__(self,N,dim):
        self.terms = {}
        '''
        Each element of inside of terms[i] is (strength,otherterm,otherlocation)
        If its a linear term, then otherterm is none and otherloc are zeros.
        '''
        self.dim = dim
        for i in range(N**2-1):
            self.terms[i] = []
        pass
    
    def add_term(self,term):
        '''
        Add a term to the Hamiltonian.
        <term> is an array with each element being a single power.
        The first element is the center index and is 2 long,
        describing the var number on the site, and the strength. All following
        indices are 1+N long, where N is the dimension; the first
        is the var number, and the other N describe the offset from
        the center.
        
        Because interactions are assumed translationally invarient,
        the symmetric counterpart is also included; eg (0,0,1) and (0,0,-1)
        
        Example term:
        
        ((0,1)) - A strength-1 local interaction on var 0
        ((0,1),(1,array(0,1))) - Strength-1 interaction in 2d between var 0 and its neighbor 1
        '''
        
        
        if len(term)==1:
            self.terms[term[0][0]].append([term[0][1],None,zeros(5)])
        elif len(term)==2:
            self.terms[term[0][0]].append([term[0][1],term[1][0],term[1][1]])
            self.terms[term[1][0]].append([term[0][1],term[0][0],-term[1][1]])
        else:
            raise 'Too uncreative to program past two. Read Comments below'
            '''
            To expand past terms beyond quadratic, one must also change dH.
            Have each bit be (J,[(term,offset),(term,offset),...]) or similar...
            '''

    
    def print_hamiltonian(self):
        '''
        Displays a graphical form of the hamiltonian.
        This is purely cosmetic...
        '''
                
        pass
    
    def dH(self,data,ind):
        '''
        Differentiates w.r.t. dH on index i and returns the relevent data.
        Returned is an array shape [ dim data ]
        '''
        dH = zeros(data.shape[0:-1])
        for term in self.terms[ind]:
            if not term[1]:
                dH +=term[0]
            else:
                #print data
                # Again, I'm doing things not in general cases... tsk tsk
                if self.dim==1:
                    dH += term[0]*roll(data[:,term[1]],term[2][0],axis=0)
                elif self.dim==2:
                    dH += term[0]*roll(roll(data[:,:,term[1]],term[2][1],axis=1),term[2][0],axis=0)
                elif self.dim==3:
                    dH += term[0]*roll(roll(roll(data[:,:,:,term[1]],term[2][2],axis=2),term[2][1],axis=1),term[2][0],axis=0)
                else:
                    raise 'Dimensionality too high! Im doing special case D<4'
        return dH


        