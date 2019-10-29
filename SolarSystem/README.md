# Solar System
Read more about this project at my [website](https://simehaa.github.io/solar-system/).
The C++ code is set up to simulate the Solar system (Sun, the eight planets and of course Pluto).
## Algorithm
The motion of the 10-body problem is solved by using the *velocity Verlet method*. Please note that armadillo (non standard C++ library) is necessary in order to run the program. If you wish to clone and download repository, feel free to do so. To install armadillo in linux, simply run
```
$ sudo apt-get install liblapack-dev libblas-dev libboost-dev libarmadillo-dev
```
## Play around with the simulations
The main.cpp file is set to initialize a solar system with initial values in the file "./data/initial_values.txt". The data that is currently there is from Nasa.gov and includes the Sun, the eight planets and Pluto with origin at mass center. Of course, any number of objects can be put there: with customizable masses, initial position and velocity. The solar system will be simulated using the total time **T** and time step **dt** in main.cpp. Construct the executable file and run it by:
```
$ make
$ ./test.x
```
Then the results can easily visualized by running:
```
$ python3 plot.py
```
In the plot.py program, it can be specified (in the main function) to plot a static plot or an animation.
