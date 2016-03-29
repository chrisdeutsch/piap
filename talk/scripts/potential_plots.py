import numpy as np
import matplotlib.pyplot as plt

plt.style.use("publication")
plt.rcParams.update({"font.size": 18})

### Coulomb with hard core
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

plt.subplots_adjust(left=0.18, bottom=0.18)
plt.savefig("figures/potential_coulomb.pdf")
plt.close()


### Lennard Jones
r = np.linspace(0.2, 3.0, 1000)

pot = 1.0 / r**12 - 1.0 / r**6

plt.plot(r, pot)

plt.grid()
plt.xlim(0.5, 2.5)
plt.ylim(-0.4, 1.0)
plt.xlabel("distance~$r$")
plt.ylabel("potential~$V$")

plt.subplots_adjust(left=0.18, bottom=0.18)
plt.savefig("figures/potential_lennard_jones.pdf")
plt.close()