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


AVNER_TIME = 3
DEFAULT_TIME = AVNER_TIME / 3E10
PATH_X = "data/y.txt"
PATH_Y = "data/y.txt"
X_MAX = 3
TH = 1160500
arad = 7.5657674E-15


def get_row_by_time(time, file):
    prev_time = 1e-14
    with open(file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                return [-1]
            else:
                data = line.split()
                if time < np.float64(data[0]):

                    return normalize_line(data)
                prev_time = data[0]


def normalize_line(data):
    if data[0] == -1:
        print ("No time found")
    data = data[:-1]
    data = np.array(data, dtype=np.float64)
    return data


def get_row_by_cycle(cycle, file):
    prev_time = 1e-14
    with open(file, "r") as f:
        lines = f.readlines()
        return lines[cycle].split(' ')


def main(path_to_data, time, cycle):
    if time == -1:
        data = get_row_by_cycle(cycle, path_to_data)
    else:
        data = get_row_by_time(time, path_to_data)
        data_r = get_row_by_time(time, "data/energy.txt")
    with open(PATH_X, "r") as f:
        x = f.readline().split()
        x = normalize_line(x)
    if data[0] == -1:
        print ("No time found")
    else:
        # because we have a \n in the end of the line we gotta truncate it
        size = min(len(x) - 1, len(data) - 2)
        # Normalize
        plt.plot(x[1:size-1], pow(data_r[2:size]/arad, 0.25) / TH, "k", label="Tr")
        plt.plot(x[1:size - 1], (data[2:size]) / TH, "--r", label="Tm")
        plt.xlim(0, X_MAX)
        plt.ylim([0, 1.5])
        xtics = np.arange(0, X_MAX, 0.5)
        ytics = np.arange(0, 1, 0.1)
        plt.legend()
        plt.xticks(xtics, xtics)
        plt.yticks(ytics, ytics)
        plt.ylabel(os.path.basename(path_to_data))
        plt.title("For time: {0}. In Avner paper: {1}".format(data[0], AVNER_TIME))
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs tests with varying input sizes.')
    parser.add_argument('-d',
                        dest='d',
                        default="data/temperature.txt",
                        help='Path to the directory containing the runs.')
    parser.add_argument('-t',
                        dest='t',
                        default=-1,
                        help='Path to the directory containing the runs.')
    parser.add_argument('-cyc',
                        dest='cyc',
                        default=-1,
                        help='Path to the directory containing the runs.')
    args = parser.parse_args()
    if args.t  == -1:
        args.t = DEFAULT_TIME
    if args.cyc != -1:
        args.t = -1
    main(os.path.abspath(args.d), np.float64(args.t), np.int(args.cyc))


