import numpy as np
import matplotlib.pyplot as plt

plt.style.use("publication")


# Coulomb with hard core
r = np.linspace(0.2, 6.0, 1000)

pot_same = 1.0 / r + 1 / r**8
pot_opp = -1.0 / r + 1 / r**8

plt.plot(r, pot_same, label="same charge")
plt.plot(r, pot_opp, label="opposite charge")


plt.grid()
plt.xlim(0.5, 6.0)
plt.ylim(-1.0, 5.0)
plt.xlabel("distance~$r$")
plt.ylabel("potential~$V$")
plt.legend(loc="best")


plt.savefig("figures/potential_coulomb.pdf")
