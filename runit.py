# -*- coding: utf-8 -*-
"""
Created on Thu Sep 01 13:19:47 2016

High-level code to run sims.

Parameters are described in a config file.
Parameters include the following. Those marked with * are not necessary.
J:          Interaction Strength
outfile:    File where data is stored
U:          Self-interaction
dim:        Dimension {1,2,3}
sies:       Size per side. <30?
IC:         Initial condition, passed directly to doIT.product_IC()
tstep:      the Timestep
T:          the Total time elapsed
tobs:       Time between observations
ncycles:    Number of times to repeat a run
*SU:        the group under; passed to hubbard_SUN(). {3,4,'old','SU3'}. If none, hubbard_SU3()
*obs:       observable. Passed to observable class. If none, 'superfluid'
*mu:        Disorder potential filename. If none, then no disorder.
*meanfield: do meanfield interaction {'t','f'}. If none, then 'f'

@author: Jonathan Wurtz
"""

from main import *


if __name__=="__main__":
    
    

    # Step 1: Load stuff from the configuration file
    if len(sys.argv)>1: # Put in place so I can run things on my laptop from F5...
        config_file = sys.argv[1]
    else:
        config_file = 'config.conf'
    
    f = open(config_file)
    line = f.readline()
    linecount = 1
    config_data = {}
    while line:
        print line
        linesplit = line.split(':')
        print linesplit
        linecount +=1
            
        config_data[linesplit[0].strip()]=linesplit[1].strip()
        
        line = f.readline()
    f.close()
    if linecount<11:
        raise 'Incomplete Config File!'
        
        
    # Check if the outfile is present... if its not, lets fail /before/ doing work.
    if 'outfile' not in config_data:
        raise 'No outfile defined!'
    
    # Output configuration data to file
    f = open(config_data['outfile'],'w')
    print >>f,'nparams:\t',len(config_data)
    for kkkk in config_data.keys():
        print >>f,kkkk+":\t",config_data[kkkk]
    f.close()
    
    
    
    
    # Step 2:Set up our smulation
    
    # Check default values...
    if 'meanfield' in config_data:
        domeanfield = (config_data['meanfield']=='t')
    else:
        domeanfield = False
    
    if 'obs' in config_data:
        obs_var = config_data['obs']
    else:
        obs_var = 'superfluid'
    
    if 'mu' not in config_data:
        config_data['mu'] = None
    
    if 'SU' not in config_data:
        params = Hubbard_SU3(int(config_data['dim']),int(config_data['sies']),double(config_data['J']),double(config_data['U']))
    else:
        if config_data['SU'].isdigit():
            params = Hubbard_SUN(int(config_data['dim']),int(config_data['sies']),double(config_data['J']),double(config_data['U']),int(config_data['SU']),domeanfield,mu_fname=config_data['mu'])
        else:
            params = Hubbard_SUN(int(config_data['dim']),int(config_data['sies']),double(config_data['J']),double(config_data['U']),config_data['SU'],domeanfield,mu_fname=config_data['mu'])
    #params['verbose']='t'
    params['obs']=observable(obs_var,int(double(config_data['T'])/double(config_data['tobs']))+3)

        
    di = doIT(params)
    
    # Step 2.5: Add more specific ICs if necessary...
    if config_data['IC'].strip()=='diffusion':
        # Hardcoded diffusion parameters...
        if params['SU']==3:
            states = [1,2]
        elif params['SU']==4:
            states = [2,3]
        else:
            raise 'Something went Wrong!'
        inICs = zeros(di.data.shape[0:-1]).astype(object)
        if int(config_data['dim'])==2:
            for i in range(int(config_data['sies'])):
                for j in range(int(config_data['sies'])):
                    if i==0 and j==0:
                        inICs[i,j] = [['z',states[0],1]]
                    else:
                        inICs[i,j] = [['z',states[1],1]]
                
        
        else:
            raise 'Dim 2 only, becase Im lazy!'
    elif config_data['IC'].strip()=='checkerboard':
        # Hardcoded diffusion parameters...
        if params['SU']==3:
            states = [1,2]
        elif params['SU']==4:
            states = [2,3]
        else:
            raise 'Something went Wrong!'
        inICs = zeros(di.data.shape[0:-1]).astype(object)
        if int(config_data['dim'])==2:
            for i in range(int(config_data['sies'])):
                for j in range(int(config_data['sies'])):
                    if round((i+j)/2.)==(i+j)/2.: # Checkerboard pattern!
                        inICs[i,j] = [['z',states[0],1]]
                    else:
                        inICs[i,j] = [['z',states[1],1]]
                        
    else:
        inICs = eval(config_data['IC'].strip())
    
    # Step 3: Run it!!
    for i in range(int(config_data['ncycles'])):
        print i
        di.obs.reset()
        di.product_IC(inICs)
        di.run(double(config_data['T']),double(config_data['tstep']),double(config_data['tobs']))
        
        # Save output to file...
        di.obs.put(config_data['outfile'])
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
else:
    print 'WTF u doin, m8?'
    raise 'run this by command line'