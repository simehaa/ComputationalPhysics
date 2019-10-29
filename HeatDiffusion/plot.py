#!/usr/bin/env python3
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':12})
from matplotlib import animation
import numpy as np
import sys


def plot_times_and_errors(params,folder):
    times = [0.1,0.2]
    dx = params[2]
    j = 0
    err_nested_list = []
    for time in times:
        fig = plt.figure()
        plt.title(r"$\Delta x = %1.2f$, time=%1.2f" % (dx,time))
        ax = plt.subplot(211)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        files = ["ex","im","cn","an"]
        labels = ["Explicit\nscheme",
                  "Implicit\nscheme",
                  "Crank-\nNicolson\nscheme",
                  "Analytic\nsolution"
                  ]
        arrays = []
        errors = []
        nx,nt,dx,dt = params
        x = np.linspace(0,nx*dx,nx+1)

        # implicit, Crank-Nicolson and analytic solution plot
        for file,l in zip(files,labels):
            if file == "an":
                data = np.fromfile("data/"+folder+file+".bin", dtype=float).reshape((2,nx+1))
                arr = data[j]
            else:
                data = np.fromfile("data/"+folder+file+".bin", dtype=float).reshape((nt+1,nx+1))
                arr = data[int(time*nt)]
            arrays.append(arr)
            plt.plot(x,arr,'.',label=l)
        ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
        plt.ylabel('Temperature [ ]')
        plt.grid()

        ax = plt.subplot(212)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        for i in range(3):
            err = np.abs(arrays[i] - arrays[3])
            plt.plot(x,err,'.')
        plt.xlabel('Position [ ]')
        plt.ylabel('Absolute error [ ]')
        plt.yscale('log')
        plt.grid()
        j = 1 # do not remove, fixes index for analytic: data file only has the two times 0.1, 0.2


if __name__=="__main__":
    # dx = 0.1
    params = [10,1000,0.1,0.001]
    plot_times_and_errors(params,"dx01/")
    # dx = 0.01
    params = [100,100000,0.01,0.00001]
    plot_times_and_errors(params,"dx001/")
    plt.show()
