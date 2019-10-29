# Numerical Simulation of Heat Diffusion
The heat diffusion is simulated in one and two dimensions in C++.
## How to run
Make sure armadillo (for Ubuntu: libarmadillo-dev, libblas-dev, liblapack-dev) is
installed. Modify src/main.cpp for what simulations you want to run, the run:
```
$ make
$ ./test.x
```
To compile and run the program. The corresponding plot functions can be found in
* plot1d.py
* plot2d.py
* plotlithosphere.py
which can be used to illustrate the results.
## Gifs: shows the time evolution temperature distribution
* twodim.gif, shows a 2D special case, which has a closed-form analytic solution.
* lithosphere_decay.gif, shows the time evolution of the lithosphere with decaying radioactive materials in the mantle.
* lithosphere_nodecay.gif, shows the time evolution of the lithosphere with radioactive materials in the mantle.
* lithosphere_decay_impact.gif, shows the time evolution the difference between lithosphere_nodecay.gif and lithosphere_decay.gif.
