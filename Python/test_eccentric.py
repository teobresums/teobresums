# Run a test
import EOBRun_module
import matplotlib.pyplot as plt
import time

#-------------------
# SXS:BBH:1355 - q=1
# omg_a = 0.02805750
#-------------------
f0   = 18.132161296614207
e0   = 0.0890000000
#-------------------
# SXS:BBH:1356 - q=1
# omg_a = 0.019077
#-------------------
f0  = 12.328512556553836
e0  = 0.1503800000
#-------------------
# SXS:BBH:1358 - q=1
# omg_a = 0.021238
#-------------------
f0  = 13.725058954557340
e0  = 0.1807800000
#-------------------
# SXS:BBH:1359 - q=1
# omg_a = 0.02139
#-------------------
f0  = 13.823288964967581
e0  = 0.18240
#-------------------
# SXS:BBH:1357 - q=1
# omg_a = 0.01960
#-------------------
f0  = 12.666501342373286
e0  = 0.19201
#-------------------
# SXS:BBH:1361 - q=1
# omg_a = 0.0210
#-------------------
f0 = 13.571251438257093
e0 = 0.23557
#-------------------
# SXS:BBH:1360 - q=1
# omg_a = 0.01959
#-------------------
f0   = 12.6600388
e0   = 0.2429
#------------------
#SXS:BBH:1362 - q=1
# omg_a = 0.01914
#------------------
f0 = 12.369226310868608
e0 = 0.3019
#------------------
#SXS:BBH:1363 - q=1
# omg_a = 0.01908
#------------------
f0 = 12.330451306759301
e0 = 0.30479
#--------------------
#SXS:BBH:1364   - q=2
# omg_a = 0.025231
#---------------------
f0 = 16.305535478031654
e0 = 0.08464
#--------------------
#SXS:BBH:1365   - q=2
# omg_a = 0.023987
#---------------------
f0 = 15.501600392832042
e0 = 0.11015
#--------------------
#SXS:BBH:1366   - q=2
# omg_a = 0.02580
#---------------------
f0 = 16.673251767001570
e0 = 0.1496
#--------------------
#SXS:BBH:1367   - q=2
# omg_a = 0.02580
#---------------------
f0 = 16.673251767001570
e0 = 0.15065
#--------------------
#SXS:BBH:1368   - q=2
# omg_a = 0.02527
#---------------------
f0 = 16.330739230702701
e0 = 0.14951
#--------------------
#SXS:BBH:1369   - q=2
# omg_a = 0.01735
#---------------------
f0 = 11.212438688274313
e0 = 0.3134
#--------------------
#SXS:BBH:1370   - q=2
# omg_a = 0.016915
#---------------------
f0 = 10.931319908481843
e0 = 0.31445
#--------------------
#SXS:BBH:1371   - q=3
# omg_a = 0.029058
#---------------------
f0 = 18.778734490136888
e0 = 0.0912
#--------------------
#SXS:BBH:1372   - q=3
# omg_a = 0.026070
#---------------------
f0 = 16.847739285493446
e0 = 0.14915
#--------------------
#SXS:BBH:1373   - q=3
# omg_a = 0.0253
#---------------------
f0 = 16.350126732757353
e0 = 0.15035
#--------------------
#SXS:BBH:1374   - q=3
# omg_a = 0.016946
#---------------------
f0 = 10.951353660604985
e0 = 0.31388


# SXS:BBH:324
# omg0 = 0.01894
# f0   = 12.239976297170919
# e0   = 0.2948000000

pars = {
    'M'                  : 50.,
    'q'                  : 3.,
    'chi1'               : 0.,
    'chi2'               : 0.,
    'Lambda1'            : 0.,
    'Lambda2'            : 0.,     
    'domain'             : 0,      #Set 1 for FD. Default = 0
    'arg_out'            : 1,      #Output hlm/hflm. Default = 0
    'use_mode_lm'        : [1],    #List of modes to use/output through EOBRunPy
    'output_lm'          : [1],    #List of modes to print on file
    'srate_interp'       : 4096.,  #srate at which to interpolate. Default = 4096.
    'use_geometric_units': 0,      #output quantities in geometric units. Default = 1
    'df'                 : 0.01,   #df for FD interpolation
    'initial_frequency'  : f0,     #in Hz if use_geometric_units = 0, else in geometric units
    'interp_uniform_grid': 1,      #interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
    'ecc'                : e0,     #Eccentricity. Default = 0.
    'ecc_freq'           : 2,      #Use periastron (0), average (1) or apastron (2) frequency for initial condition computation. Default = 1
}

#Run the WF generator
start = time.time()
t, hp, hcm, hlm, dyn = EOBRun_module.EOBRunPy(pars)
end = time.time()
DeltaT = end-start
print("Full time=%s"%DeltaT);

plt.plot(t, hp)
plt.show()
