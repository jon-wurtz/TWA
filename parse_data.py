# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 13:34:22 2016

@author: Jonathan Wurtz
"""

# Parses output files which are copied here...
import os
from numpy import *


filename = 'small_mottIC.o57511-'
os.chdir('C:\Users\Jonathan Wurtz\Documents\Research\Scripts\SU3_dynamics\Outputs_0826')
filed = os.listdir(os.getcwd())

f = open(filename+'0')



header = []
for i in range(8):
    header.append(f.readline())
    
ff = f.readline()
datstring = ''
while ff:
    if 'D' not in ff and '[' not in ff:
        datstring +=ff
    ff = f.readline()
    
dat = genfromtxt(StringIO(datstring))
f.close()
n_tsteps = 1.0*dat.shape[0]/int(header[6])

if round(n_tsteps)!=n_tsteps:
    raise 'Bad!'
    
dat_out = dat.reshape(int(header[6]),n_tsteps,3).sum(0)/int(header[6])
kk = 1


for k in filed:
    if filename in k:
        print 'ping'
        f = open(k)
        for i in range(8):
            f.readline()
        
        ff = f.readline()
        datstring = ''
        while ff:
            if 'D' not in ff and '[' not in ff:
                datstring +=ff+str('\n')
            ff = f.readline()
        dat_out += genfromtxt(StringIO(datstring)).reshape(int(header[6]),n_tsteps,3).sum(0)/int(header[6])
        kk+=1
dat_out = dat_out/kk
figure()
plot(dat_out[:,0],dat_out[:,2])
title('Bose Hubbard Model, IC: '+header[5]+'(J,U)=('+header[1][0:-1]+', '+header[2][0:-1]+')\n Number of Trials:'+repr(kk*int(header[6]))+'\n tstep: '+header[7][0:-1]+'\nData:'+filename)
xlabel('Scaled Time tU')
ylabel('Order Parameter')
