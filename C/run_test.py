# Run a test
import EOBRunTD_module
import matplotlib.pyplot as plt

t, hp, hc = EOBRunTD_module.EOBRunTD(2.7, 1., 0., 0., 100.,100.)

plt.plot(t, hp)
plt.show()
