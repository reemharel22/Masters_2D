import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax
import matplotlib.animation as animation
import time
import argparse
import os
import re
import subprocess
import distutils
from threading import Thread
import time
import threading
import shutil
DEFAULT_TIME = 1E-2
PATH_X = "data/x.txt"
X_MAX = 2


def get_row_by_time(time, file):
    prev_time = 1e-14
    with open(file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                return [-1]
            else:
                data = line.split(' ')
                if time < np.float64(data[0]):
                    return data
                prev_time = data[0]


def main(path_to_data, time):
    with open(PATH_X, "r") as f:
        x = f.readline().split()
    data = get_row_by_time(time, path_to_data)
    if data[0] == -1:
        print ("No time found")
    else:
        x = x[:-1]
        data = data[:-1]
        data = np.array(data, dtype=np.float64)
        x = np.array(x, dtype=np.float64)
        # because we have a \n in the end of the line we gotta truncate it
        size = min(len(x) - 1, len(data) - 2)
        plt.plot(x[:size-1], data[1:size], "r")
        plt.xlim(0, X_MAX)
        plt.ylim([0, 1])
        xtics = np.arange(0, X_MAX, 0.5)
        ytics = np.arange(0, 1, 0.1)
        plt.xticks(xtics,xtics)
        plt.yticks(ytics,ytics)
        plt.ylabel(os.path.basename(path_to_data))
        plt.title("For time: {0} ".format(data[0]))
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs tests with varying input sizes.')
    parser.add_argument('-d',
                        dest='d',
                        default="data/energy.txt",
                        help='Path to the directory containing the runs.')
    parser.add_argument('-t',
                        dest='t',
                        default=-1,
                        help='Path to the directory containing the runs.')
    args = parser.parse_args()
    if args.t  == -1:
        args.t = DEFAULT_TIME
    main(os.path.abspath(args.d), np.float64(args.t))


