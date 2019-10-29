#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams.update({"font.size": 18})
import numpy as np
import sys


def plot():
    """
    Output:
    Plots the four numerical results and three analytic results.
    Plots the error between the three cases with both analytic and numerical
    results.
    """
    files = [
        "data/noheat/u100.bin",
        "data/naturalheat/u100.bin",
        "data/nodecay/u100.bin",
        "data/decay/u100.bin",
    ]
    labels = [
        "Numerical: (1) no heat",
        "Numerical: (2) natural heat",
        "Numerical: (3) enriched mantle, no decay",
        "Numerical: (4) enriched mantle with decay",
    ]
    colors = ["b-", "y-", "r-", "g-"]
    numerical = []
    nx = 126
    ny = 101
    nt = 100
    fig = plt.figure()
    ax = plt.subplot(111)
    for file, label, c in zip(files, labels, colors):
        image = np.fromfile(file, dtype=float).reshape((ny, nx))
        temp = image[:, int(nx / 2)] * 1292 + 8
        numerical.append(temp)
        depth = np.linspace(0, 120, ny)
        plt.plot(depth, temp, c, label=label)
    plt.xlabel("Depth [km]")
    plt.ylabel(r"Temperature $[^\circ C]$")
    plt.grid()

    # analytic: no heat production
    analytic = []
    x = np.linspace(0, 120, 101)
    y = np.linspace(0, 1, 101) * 1292 + 8
    analytic.append(y)
    plt.plot(x, y, "b--", label="Analytic: (1) no heat")

    # analytic: natural heat production
    xzones = [x[:17], x[17:34], x[34:]]
    zones = [[-0.28, -23.66, 8], [-0.07, -15.26, 92], [-0.01, -10.46, 188]]
    y2 = []
    for x, zone in zip(xzones, zones):
        y = np.polyval(zone, -x)
        y2.append(y)
        plt.plot(x, y, "y--", label="Analytic: (2) natural heat")
    analytic.append(np.concatenate((y2[0], y2[1], y2[2])))
    # analytic: natural + subduction
    zones = [[-0.28, -29, 8], [-0.07, -20.6, 92], [-0.11, -23.8, 28]]
    y3 = []
    for x, zone in zip(xzones, zones):
        y = np.polyval(zone, -x)
        y3.append(y)
        plt.plot(x, y, "r--", label="Analytic: (3) enriched mantle, no decay")
    analytic.append(np.concatenate((y3[0], y3[1], y3[2])))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + 0.3, box.width, box.height * 0.7])
    ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, -0.75), fancybox=True, shadow=True
    )
    plt.show()

    # error plot
    labels = ["(1) no heat", "(2) natural heat", "(3) enriched mantle,\n      no decay"]
    fig = plt.figure()
    ax = plt.subplot(111)
    x = np.linspace(0, 120, 101)
    for a, n, lab in zip(analytic, numerical, labels):
        plt.plot(x[1:-1], np.abs(a[1:-1] - n[1:-1]) / a[1:-1], label=lab)
    plt.yscale("log")
    plt.xlabel("Depth [km]")
    plt.ylabel("Relative error")
    plt.grid()
    ax.legend(loc=3, fancybox=True, shadow=True)
    plt.show()


def animate():
    """
    Animates and saves two gifs: temperature evolution for the simulations
    with enriched mantle with and without decay
    """
    infiles = ["data/nodecay/u", "data/decay/u"]
    outfiles = ["lithosphere_nodecay.gif", "lithosphere_decay.gif"]
    nx = 126
    ny = 101
    nt = 100
    for infile, outfile in zip(infiles, outfiles):
        images = []
        for i in range(1, nt):
            image = np.fromfile(infile + str(i) + ".bin", dtype=float)
            images.append(image.reshape((ny, nx)) * 1292 + 8)

        fig = plt.figure()
        im = plt.imshow(images[0], animated=True)

        def updatefig(j):
            im.set_array(images[j])
            return [im]

        ani = animation.FuncAnimation(
            fig, updatefig, frames=range(nt - 1), interval=70, blit=True
        )

        plt.colorbar()
        plt.xlabel(r"$n_x$")
        plt.ylabel(r"$n_y$")
        ani.save(outfile, writer="imagemagick", fps=30)


def plot_difference_nodecay_decay():
    """
    Plots the difference at the end of the simulation with and without decay
    """
    fileDecay = "data/decay/u100.bin"
    fileNoDecay = "data/nodecay/u100.bin"
    nx = 126
    ny = 101
    nt = 100
    differenceImages = []
    plt.rcParams.update({"font.size": 15})
    decay = np.fromfile(fileDecay, dtype=float).reshape((ny, nx)) * 1292 + 8
    noDecay = np.fromfile(fileNoDecay, dtype=float).reshape((ny, nx)) * 1292 + 8
    diff = noDecay - decay
    x = np.linspace(0, 150, nx)
    y = np.linspace(0, 120, ny)
    levels = np.linspace(0, 35, 100)
    plt.contourf(x, y, diff, levels=levels)
    cbar = plt.colorbar(ticks=[0, 7, 14, 21, 28, 35])
    cbar.set_label(r"$T_{no decay}(t) - T_{decay}(t), [^\circ C]$")
    plt.xlabel("Width [km]")
    plt.ylabel("Depth [km]")
    plt.show()

    maxtemp = np.max(diff)
    depth = np.argmax(diff[:, int(nx / 2)], axis=0)
    print("max temperature difference: %f, depth: %f" % (maxtemp, depth * 1.2))


def animate_difference_nodecay_decay():
    """
    Animate the time evolution of the temperature difference for the two
    simulations with and without decay.
    """
    fileDecay = "data/decay/u"
    fileNoDecay = "data/nodecay/u"
    nx = 126
    ny = 101
    nt = 100
    differenceImages = []
    plt.rcParams.update({"font.size": 10})
    for i in range(1, nt):
        decay = (
            np.fromfile(fileDecay + str(i) + ".bin", dtype=float).reshape((ny, nx))
            * 1292
            + 8
        )
        noDecay = (
            np.fromfile(fileNoDecay + str(i) + ".bin", dtype=float).reshape((ny, nx))
            * 1292
            + 8
        )
        differenceImages.append(noDecay - decay)

    fig = plt.figure()
    im = plt.imshow(differenceImages[0], vmin=0, vmax=40.0, animated=True)

    def updatefig(j):
        im.set_array(differenceImages[j])
        return [im]

    ani = animation.FuncAnimation(
        fig, updatefig, frames=range(nt - 1), interval=50, blit=True
    )

    cbar = plt.colorbar()
    cbar.set_label(r"$T_{no decay}(t) - T_{decay}(t), [^\circ C]$")
    plt.title("Decay impact\non temperature")
    plt.xlabel(r"$n_x$")
    plt.ylabel(r"$n_y$")
    ani.save("lithosphere_decay_impact.gif", writer="imagemagick", fps=30)


if __name__ == "__main__":
    plot()
    animate()
    plot_difference_nodecay_decay()
    animate_difference_nodecay_decay()
