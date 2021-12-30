import pylab
import pycbc.waveform
from pycbc.waveform import td_approximants

# This tells pycbc about our new waveform so we can call it from standard
# pycbc functions. If this were a frequency-domain model, select 'frequency'
# instead of 'time' to this function call.
#
#from teobresums import teobresums_td
#pycbc.waveform.add_custom_waveform('teobresums', teobresums_td, 'time', force=True)

print(td_approximants())
 
# Let's plot what our new waveform looks like
hp, hc = pycbc.waveform.get_td_waveform(approximant="teobresums",
                                        mass1=1.4, mass2=1.4, lambda1=400.,lambda2=400., spin1z=0., spin2z=0.,
                                        f_lower=20.,
                                        delta_t=1/2048.)#delta_f=1./128)
pylab.figure(0)
pylab.plot(hp.sample_times, hp)
pylab.xlabel('Time (s)')

pylab.figure(1)
hf = hp.to_frequencyseries()
pylab.plot(hf.sample_frequencies, hf.real())
pylab.xlabel('Frequency (Hz)')
#pylab.xscale('log')
#pylab.xlim(20, 100)
pylab.show()

