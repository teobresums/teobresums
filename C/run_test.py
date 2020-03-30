# Run a test
import EOBRun_module
import matplotlib.pyplot as plt

pars = {
    'M': 50.,
    'q': 1.,
    'chi1': 0.,
    'chi2': 0.,
    'Lambda1': 0.,
    'Lambda2': 0.,
    'interp_uniform_grid':0,
    'use_geometric_units':0,
    'initial_frequency':20,
    'ecc':0.1
    }

t, hp, hc = EOBRun_module.EOBRunTD(pars)

plt.plot(t, hp)
plt.show()
