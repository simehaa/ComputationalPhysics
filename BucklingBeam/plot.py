#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import re

color = ["r--", "g--", "b--", "m--"]
for i in range(4):
    rho = []
    psi = []
    with open("./data/project2_" + str(i) + ".txt", "r") as infile:
        # first three lines
        line = infile.readline()
        words = re.sub(
            ":,[^\w]", " ", line
        ).split()  # splits line into list of words that were seperated by : or , or character
        eigval = eval(words[1])

        line = infile.readline()
        words = re.sub(":,[^\w]", " ", line).split()
        wr = eval(words[1])

        infile.readline()

        # the rest of the lines (rho_i psi_i):
        for line in infile:
            r, p = line.split()
            rho.append(eval(r))
            psi.append(eval(p))
        infile.close()

    psi = np.asarray(psi)
    psi *= psi
    lab = r"$\lambda = %5.2f, \omega_r = %4.2f$" % (eigval, wr)
    plt.plot(rho, psi, color[i], label=lab)

plt.title("Two electron in a harmonic oscillator with coulomb interaction")
plt.xlabel(r"natural length scale: $\rho$")
plt.ylabel(r"probability distribution: $|\psi|^2$")
plt.grid()
plt.legend()
plt.savefig("./data/project2.png")
