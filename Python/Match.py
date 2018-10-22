from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

LAL_MTSUN_SI    = 4.925491025543575903411922162094833998e-6
verbose         = 1

# Total masses in units of M_SUN 
M_tot  = [10,60]
# Data to form the file name
q      = [1,5]
freq   = ['0001', '0006']
spin1  = ['00001','01','01','-01','099','099','-099']
spin2  = ['00001','01','-01','-01','099','-099','-099']

# psd
f, psd  = np.loadtxt("/mnt/c/Users/giuli/src/lalsuite/lalsimulation/src/LIGO-P1200087-v18-AdV_DESIGN.txt",unpack=True) #asd
psd    *= psd
# generate an interpolant for the PSD
psd_int =  interp1d(f, psd, bounds_error=False, fill_value=np.inf)


match  = []
for i in range(len(q)):

    for j in range( len(spin1) ):
        
        #hc has a minus sign with respect to the one unpacked here
        t_eob, hp_eob, hc_eob = np.loadtxt("/mnt/c/Users/giuli/Repositories/teobresums/C/bbh_q{0}_s{1}s{2}_f{3}/waveform.txt".format(q[i], spin1[j], spin2[j], freq[i]), unpack=True)
        t_pa, hp_pa, hc_pa    = np.loadtxt("/mnt/c/Users/giuli/Repositories/teobresums/C/bbh_q{0}_s{1}s{2}_f{3}_postadiab/waveform.txt".format(q[i], spin1[j], spin2[j], freq[i]), unpack=True)

        M               = M_tot[i]*LAL_MTSUN_SI
         
        dt              =  0.5*M                #from parfile
        n1              = len(hp_eob)
        n2              = len(hp_pa)
        n               = np.minimum(n1,n2) 


        # Fourier transforms
        hptilde_eob  = np.fft.rfft(hp_eob[:n])*dt# we multiply by dt to get the dimensionfull FFT
        hctilde_eob  = np.fft.rfft(-hc_eob[:n])*dt# real hc has a minus sign

        hptilde_pa   = np.fft.rfft(hp_pa[:n])*dt
        hctilde_pa   = np.fft.rfft(-hc_pa[:n])*dt
        # Time shifts and normalisation of the Fourier transform                       
        tc_index_eob    = np.argmax(hp_eob[:n]**2+hc_eob[:n]**2)       
        tc_eob          = t_eob[tc_index_eob]
        tc_index_pa     = np.argmax(hp_pa[:n]**2+hc_pa[:n]**2)       
        tc_pa           = t_eob[tc_index_pa]
        timeShift       =(tc_pa - tc_eob)*M # Arbitrary time shift
        frequencies     = np.fft.rfftfreq(n,d=dt)/M
        sigmasq         =  psd_int(frequencies) * dt**2
        # Antenna pattern        
        fp=1
        fc=0    
        # Projection on the detector
        template      = (fp*hptilde_pa+fc*hctilde_pa)*np.exp(-1j*2.0*np.pi*timeShift*frequencies)
        data          = (fp*hptilde_eob+fc*hctilde_eob)
        # Match
        overlap       = np.sum( (np.conj(data)*template+np.conj(template)*data)/sigmasq ).real
        norm1         = np.sum( 2*np.conj(template)*template/sigmasq ).real
        norm2         = np.sum( 2*np.conj(data)*data/sigmasq ).real
        Normalisation = np.sqrt(  norm1*norm2 )
        temp_match    = overlap/Normalisation
        
        match.append(temp_match)
        
        if (verbose):
            print temp_match

      