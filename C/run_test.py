# Run a test
import EOBRun_module
import matplotlib.pyplot as plt

# Example un

pars = {
    'M': 50.,
    'q': 1.2,
    'chi1': 0.,
    'chi2': 0.,
    'Lambda1': 0.,
    'Lambda2': 0.,
    'domain':0,      #Set 1 for FD
}

t, hp, hcm, hlm = EOBRun_module.EOBRunPy(pars)
#f, hpr, hpi, hcr, hci, hflm = EOBRun_module.EOBRunPy(pars) # If FD

# Plots
plt.plot(t, hp)
plt.show()
plt.plot(t, hlm['1'][0]) #22 Amplitude
plt.show()
