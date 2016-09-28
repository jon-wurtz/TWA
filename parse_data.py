# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 13:34:22 2016

@author: Jonathan Wurtz
"""

# Parses output files which are copied here...
import os
from numpy import *

def new_parse2(filename,header_only=False):

    os.chdir('C:\Users\Jonathan Wurtz\Documents\Research\Scripts\SU3_dynamics\Outputs_0826')
    filed = os.listdir(os.getcwd())
    touse = []
    # Find number of files of the same type...
    for i in range(256): # Ug, lazy solution... what if you ahve more then 256 jobs?
        if filename.replace('*',str(i)) in filed and filename.replace('*',str(i)) not in touse:
            touse.append(filename.replace('*',str(i)))
        
        
    f = open(touse[0],'rb')
    
    data = f.read().replace('\r\n','\n').split('---DATA---')
    
    config_data = {}
    for line in data[0].split('\n')[0:-1]:
        config_data[line.split(':\t')[0]] = line.split(':\t')[1]
    
    if header_only==True:
        for k in config_data:
            print k,str('\t'),config_data[k]
        return None,None,config_data
        
            
    dat_out_list = []
    for dats in data[1::]:
        try:    
            dat_out_list.append(fromstring(dats[1:-1]))
        except:
            pass
    
    f.close()
    
    for touse_ in touse[1::]:
        #print touse_
        f = open(touse_,'rb')
        data = f.read().replace('\r\n','\n').split('---DATA---')
        f.close()
        for dats in data[1::]:
            try:
                dat_out_list.append(fromstring(dats[1:-1]))
            except:
                pass
    
    dat_out_arr = array(dat_out_list)
    config_data['nruns'] = len(dat_out_list)
    
    num_obs = int(double(config_data['T'])/double(config_data['tobs']))+1
    
    
    if len(dat_out_list)==0:
        print 'No data!'
        return None,None,config_data
    TT = linspace(0,double(config_data['T']),len(dat_out_list[0]))
    if abs(len(dat_out_list[0])-num_obs)<3: #1d data
        #print average(dat_out_arr,0).shape
        #print TT.shape
        plot(TT,average(dat_out_arr,0),'b',linewidth=2)
        plot(TT,average(dat_out_arr,0)+std(dat_out_arr,0)/sqrt(len(dat_out_list)),'r--')
        plot(TT,average(dat_out_arr,0)-std(dat_out_arr,0)/sqrt(len(dat_out_list)),'r--')
    elif int(sqrt(len(dat_out_list[0])))==num_obs: #2d data
        pass
    
    
        
    return dat_out_arr,TT,config_data

def scan_files(filename,fnumbers,keys=None):
    # Filename replaced by " # "
    os.chdir('C:\Users\Jonathan Wurtz\Documents\Research\Scripts\SU3_dynamics\Outputs_0826')
    filename = filename.replace('*','0') # Pick the zeroth file
    
    filed = os.listdir(os.getcwd())
    
    for fnum in fnumbers:
        if filename.replace('#',str(fnum)) in filed:
            f = open(filename.replace('#',str(fnum)),'rb')
            data = f.read().replace('\r\n','\n').split('---DATA---')
            f.close()    
            
            config_data = {}
            for line in data[0].split('\n')[0:-1]:
                config_data[line.split(':\t')[0]] = line.split(':\t')[1]
                
            print fnum
            if keys==None:
                for k in config_data:
                    print k,str('\t'),config_data[k]
            else:
                for k in keys:
                    if k in config_data:
                        print k,str('\t'),config_data[k]
            
'''

#return dat_temp    


    
    
    
'''