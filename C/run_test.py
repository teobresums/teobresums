# Run a test
import EOBRunTD_module
import matplotlib.pyplot as plt

pars = {
    'M': 50.,
    'q': 1.,
    'chi1': 0.,
    'chi2': 0.,
    'Lambda1': 0.,
    'Lambda2': 0.
}

t, hp, hc = EOBRunTD_module.EOBRunTD(pars)

plt.plot(t, hp)
plt.show()
