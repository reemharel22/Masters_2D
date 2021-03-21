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
PATH_2 = "data/x.txt"
PATH_ENERGY = "data/energy.txt"
PATH_TEMP = "data/temperature.txt"
PATH_ENERGY1 = "data/energy1.txt"
PATH_TEMP1 = "data/temperature1.txt"
X_MAX = 3
TH = 1160500
arad = 7.5657674E-15
COLORS=["r", "k", "b"]


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
    t = data[0]
    data = data[:-1]
    data = np.array(data, dtype=np.float64)
    return data


def get_row_by_cycle(cycle, file):
    prev_time = 1e-14
    with open(file, "r") as f:
        lines = f.readlines()
        return lines[cycle].split(' ')


def prepare_data_for_plot(time, two_plots=False):
    data_temp = []
    data_energy = []

    data_temp1 = get_row_by_time(time, PATH_TEMP)
    data_energy1 = get_row_by_time(time, PATH_ENERGY)
    data_temp.append(data_temp1)
    data_energy.append(data_energy1)
    if two_plots:
        data_temp2 = get_row_by_time(time, PATH_TEMP1)
        data_energy2 = get_row_by_time(time, PATH_ENERGY1)
        data_temp.append(data_temp2)
        data_energy.append(data_energy2)

    for i in range(len(data_temp)):
        for j in range(len(data_temp[i])):
            if j != 0:
                data_temp[i][j] = (float(data_temp[i][j]) / TH)
                data_energy[i][j] = pow(data_energy[i][j] / arad, 0.25) / TH
    with open(PATH_X, "r") as f:
        x = f.readline().split()
        x = normalize_line(x)
    # because we have a \n in the end of the line we gotta truncate it
    size = min(len(x) - 1, len(x) - 2)
    # Normalize
    x = x[1:size-1]
    t = data_energy[0][0]
    data_temp_final = []
    data_energy_final = []
    for i in range(len(data_temp)):
        data_temp_final.append(data_temp[i][2:size])
        data_energy_final.append(data_energy[i][2:size])
    return x, t, data_temp_final, data_energy_final


def plot_surf():
    

def main(time, num_second_plot):
    if num_second_plot == -1:
        two_plots = False
        legend_lbl = [""]
    else:
        two_plots = True
        with open(PATH_2, "r") as f:
            lines = f.readline()
            print lines
            legend_lbl = ["", lines[1 + int(num_second_plot)]]
    x, t, data_temp, data_energy = prepare_data_for_plot(time, two_plots)
    for i in range(len(data_temp)):
    # Normalize
        plt.plot(x, data_energy[i][:], COLORS[i], label="Tr" + legend_lbl[i])
        plt.plot(x, data_temp[i][:], "--" + COLORS[i], label="Tm" + legend_lbl[i])
    plt.xlim(0, X_MAX)
    plt.ylim([0, 1.5])
    xtics = np.arange(0, X_MAX, 0.5)
    ytics = np.arange(0, 1, 0.1)
    plt.legend()
    plt.xticks(xtics, xtics)
    plt.yticks(ytics, ytics)
    plt.ylabel("VALUE")
    plt.title("For time: {0}. In Avner paper: {1}".format(t, AVNER_TIME))
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
    parser.add_argument('-plots',
                        dest='plots',
                        default=-1,
                        help='Path to the directory containing the runs.')
    args = parser.parse_args()
    if args.t  == -1:
        args.t = DEFAULT_TIME
    if args.cyc != -1:
        args.t = -1
    main(np.float64(args.t), args.plots)


