# -*- coding: utf-8 -*-
"""
Created on Thu Sep 01 13:19:47 2016

High-level code to run sims

@author: Jonathan Wurtz
"""

from main import *


if __name__=="__main__":
    
    # Step 1: Load stuff from the configuration file
    config_file = sys.argv[1]
    
    
    f = open(config_file)
    line = f.readline()
    linecount = 1
    config_data = {}
    while line:
        linesplit = line.split(':')
        linecount +=1
            
        config_data[linesplit[0].strip()]=linesplit[1].strip()
        
        line = f.readline()
    f.close()
    if linecount<11:
        raise 'Incomplete Config File!'
        
        
    # Check if the outfile is present... if its not, lets fail /before/ doing work.
    if 'outfile' not in config_data:
        raise 'No outfile defined!'
    '''
    print 'J:',J
    print 'U:',U
    print 'dim:',dim
    print 'sies:',sies
    print 'IC:',IC
    print 'tstep:',tstep
    print 'T:',T
    print 'tobs:',tobs
    print 'ncycles:',ncycles
    print 'outfile:',outfile
    '''
    
    # Output configuration data to file
    f = open(config_data['outfile'],'w')
    print >>f,'nparams:\t',len(config_data)
    for kkkk in config_data.keys():
        print >>f,kkkk+":\t",config_data[kkkk]
    f.close()
    
    
    
    
    # Step 2:Set up our smulation
    params = Hubbard_SU3(int(config_data['dim']),int(config_data['sies']),double(config_data['J']),double(config_data['U']))
    params['verbose']='f'
    params['obs']=observable('superfluid',int(double(config_data['T'])/double(config_data['tobs']))+1)

        
    di = doIT(params)
    
    
    # Step 3: Run it!!
    for i in range(int(config_data['ncycles'])):
        print i
        di.obs.reset()
        di.product_IC(eval(config_data['IC'].strip()))
        di.run(double(config_data['T']),double(config_data['tstep']),double(config_data['tobs']))
        
        # Save output to file...
        f = open(config_data['outfile'],'a')
        for j in range(len(di.obs.data)):
            print >>f,di.obs.T[j],di.obs.data[j]
        print >>f,'--- End of Run ---'
        f.close()
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
else:
    print 'WTF u doin, m8?'
    raise 'run this by command line'