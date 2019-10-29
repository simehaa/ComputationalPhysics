#!/usr/bin/python3
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

plt.rcParams.update({"font.size": 12})
import numpy as np


def readfile(file_number):
    """read data from one data file,
    save, and return name of planet and
    arrays of x, y, z and t"""
    filename = "./data/planet" + str(file_number) + ".txt"
    xl = []
    yl = []
    zl = []
    tl = []
    with open(filename, "r") as infile:
        planet = infile.readline().split()[0]
        print(planet)
        infile.readline()
        for line in infile:
            xval, yval, zval, tval = line.split()
            xl.append(eval(xval))
            yl.append(eval(yval))
            zl.append(eval(zval))
            tl.append(eval(tval))
        infile.close()
    x = np.asarray(xl)
    y = np.asarray(yl)
    z = np.asarray(zl)
    t = np.asarray(tl)
    return planet, x, y, z, t


def planetdata(number_of_planets):
    """generate an array data, set up
    correctly for a 3d animation"""
    data = []
    names = []
    for i in range(number_of_planets):
        planet, x, y, z, t = readfile(i)
        data_one_planet = [x, y, z]
        data.append(data_one_planet)
        names.append(planet)
    return np.asarray(data), names


def plot_3D(data, names):
    """plot a 3D plot of the number of planets
    specify limit (in AU) in the limit list"""
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    limit = [-35, 35]
    ax.set_xlim3d(limit)
    ax.set_ylim3d(limit)
    ax.set_zlim3d(limit)
    ax.set_xlabel("X [AU]")
    ax.set_ylabel("Y [AU]")
    ax.set_zlabel("Z [AU]")
    for i in range(number_of_planets):
        x, y, z = data[i]
        ax.plot(x, y, z, label=names[i])
    plt.legend()
    plt.show()


def animate_3D(data, names):
    """3D animation of all planets"""
    N = len(data[0][0])
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    limit = [-35, 35]  # AU
    ax.set_xlim3d(limit)
    ax.set_ylim3d(limit)
    ax.set_zlim3d(limit)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    def update_lines(num, dataLines, lines):
        for line, data in zip(lines, dataLines):
            line.set_data(data[0:2, :num])
            line.set_3d_properties(data[2, :num])
            # line.set_marker("o") # if not used: orbits are represented as lines
        return lines

    lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]
    ani = animation.FuncAnimation(
        fig, update_lines, N, fargs=(data, lines), interval=1, blit=False
    )

    ani.save("solar_system.gif", writer="imagemagick", fps=30)


def main():
    """set number of planet that wish to be plotted.
    data is obtained from './data/planet*.txt', where
    * = 0,1,2,..., number_of_planets - 1. Choose to plot
    the complete orbits of all planets or animate."""
    data, names = planetdata(number_of_planets)
    plot_3D(data, names)
    # animate_3D(data[::24],names)


if __name__ == "__main__":
    number_of_planets = 10
    main()

"""Terminal run:
simen@simen-ubuntu:~/steinngithub/FYS3150/Project3/OO_solar_system$ python3 plot.py
Sun
Mercury
Venus
Earth
Mars
Jupiter
Saturn
Uranus
Neptune
Pluto
"""
