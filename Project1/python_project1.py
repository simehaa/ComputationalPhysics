#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys

exponent = int(sys.argv[1])

maxerrors = []
h_list = []

for i in range(1, exponent+1):
	n = int(10**i)
	h = -i
	infile = open("results%i.txt" % i, "r")
	steps = [0]
	exact = [0]
	comp  = [0]
	err   = []
	for line in infile:
		a, b, c, d = line.split()
		steps.append(eval(a))
		exact.append(eval(b))
		comp.append(eval(c))
		err.append(eval(d))
	steps.append(1)
	exact.append(0)
	comp.append(0)

	maxerrors.append(max(err))
	h_list.append(h)

	"""
	plt.plot(steps, exact, "r.", label="Analytical solution")
	plt.plot(steps, comp, "b.", label="Numerical solution")
	plt.title("LU-decomposition solution for $n=10^%i$\n compared to the analytical expression." %i)
	plt.legend(loc="best")
	plt.ylabel("f(x) [$C$]")
	plt.xlabel("x [$m$]")
	plt.grid()
	plt.show()
	"""


plt.plot(h_list, maxerrors, "ro")
plt.title("Logarithmic scale of max error as a function of the step length.")
plt.xlabel("$log_{10}(h)$")
plt.ylabel("$log_{10}(\\epsilon)$")
plt.grid()
plt.show()
