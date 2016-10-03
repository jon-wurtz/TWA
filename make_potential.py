# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 16:44:48 2016

@author: Jonathan Wurtz
"""

import sys
from numpy import *

    




def make_random_potential(dim,sies,npix_corrlength,strength,fname,plotit=False):
    L = 10./(npix_corrlength+1e-10)
    dat = random.normal(size=(sies*ones(dim)))
    
    if dim==2:
        xx = meshgrid(linspace(-10,10,sies),linspace(-10,10,sies))
        mask = exp(-(L**-2)*(xx[0]**2 + xx[1]**2))
    if dim==3:
        xx = meshgrid(linspace(-10,10,sies),linspace(-10,10,sies),linspace(-10,10,sies))
        mask = exp(-(L**-2)*(xx[0]**2 + xx[1]**2 + xx[2]**2))
    if dim==1:
        xx = meshgrid(linspace(-10,10,sies))
        mask = exp(-sqrt((L**-2)*(xx[0]**2)))
    
    datF = real(fft.ifftn(fft.fftshift(mask)*fft.fftn(dat)))  
    # Normalize...
    datF = strength*datF/std(datF)
    
    if dim==2 and plotit:
        imshow(datF,interpolation='none')
    if dim==1 and plotit:
        #plot(datF)
        corr = zeros(int(sies/2))
        for i in range(int(sies/2)):
            corr[i] = dot(datF,roll(datF,shift=i,axis=0))
        plot(corr)
    
    f = open(fname,'w')
    savetxt(f,datF) # Retrieve with loadtxt(f) and f = open(fname,'rb')
    f.close()
    return datF
    

# sies,dim,corrlen,strength,filename
if __name__=="__main__":
    if len(sys.argv)==1:
        sies = 20               # Number of elements to a side
        dim = 2                 # Dimension; in {1,2,3}
        npix_corrlength = 2     # Correlation length (gaussian correlated)
        strength = 1            # Size of RMS fluctuations
        fname = 'test.dat'
        print 'using defaults!'
        datF = make_random_potential(dim,sies,npix_corrlength,strength,fname)
    else:
        sies = int(sys.argv[1])
        dim = int(sys.argv[2])
        npix_corrlength = double(sys.argv[3])
        strength = double(sys.argv[4])
        fname = sys.argv[5]
        
        datF = make_random_potential(dim,sies,npix_corrlength,strength,fname)    
        