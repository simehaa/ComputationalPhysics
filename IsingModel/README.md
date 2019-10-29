# Simulation of the Ising Model
This is a project in the course FYS3150, Computational Physics at University
of Oslo. The project is a numerical simulation of the Ising model and is written
by Steinn Hauser and Simen Haapnes.

For reproducability, edit what parts of the code you want to run
by commenting/uncommenting parts of src/main.cpp. All data is included in the
data folder. For edit of plots simply edit the file plot.py.
## Notes
In order to compile the project, MPI is required (download and install it first).
The program can be comiled by the simple command:
```
$ make
```
And run by
```
$ mpirun -n 4 ./test.x
```
where the number 4 represent the number of cores, and it can be switched out by any possible
number of cpu cores your computer has. If you wish to run any other part of main.cpp that is
not the study of phase transitions, it can be run by
```
$ ./test.x
```
plot.py can by run by
```
$ python3 plot.py
```
