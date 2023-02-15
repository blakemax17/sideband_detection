import pandas as pd
import numpy as np
import os
import time
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from mycolorpy import colorlist as mcp
from matplotlib.colors import ListedColormap
from scipy.optimize import curve_fit
from matplotlib.collections import LineCollection
dirname = os.path.dirname(__file__)

tic = time.perf_counter()
print(tic)

datafile = open("dmtracks.csv", "r")
df = pd.read_csv(datafile, sep=',')
df["Mainband"] = "Unclassified"
print(df["Mainband"].head(20))

events = df.groupby("initial_field")
print(events.ngroups, " fields identified")
# print(events["Mainband"].head())
# events.to_csv("Random")


def abline(slope, x0, y0, y_distance):
    """Assemble 2 points both above an below a line
    that lie along its normal vector"""
    y_above = []
    x_above = []
    for i in range(-2, 0):
        y = y0 + (i * y_distance)
        y_above.append(y0 + (i * y_distance))
        x = x0 + ((y0 - y) / slope)
        x_above.append(x)

    y_below = []
    x_below = []
    for i in range(1, 3):
        y = y0 + (i * y_distance)
        y_below.append(y0 + (i * y_distance))
        x = x0 + ((y0 - y) / slope)
        x_below.append(x)

    return [x_above, y_above, x_below, y_below]


def distance(a, b):
    """Pythagorean distance of 2 points"""
    return np.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


def is_between(a, c, b):
    """Checks if the distance from a to c and c to b
    sum to the distance from a to b"""
    accb_distance = distance(a, c) + distance(c, b)
    # print(accb_distance - distance(a, b))
    return (accb_distance >= (0.9 * distance(a, b))) & (accb_distance <= (1.1 * distance(a, b)))


iteration = 0
# print('Length:', len(df))

# Adjacent event detection
# print(events.head(20))

# Checks each track at its beginning, midoint, and endpoint
# for adjacent tracks 1 or 2 axial distances away and classifies
# the track as a mainband, dieband, or as an error
for j, field in events:
    # print(j)
    # print(field)
    # print(type(field))
    for index, event in field.iterrows():
        # print(index)
        line1 = (event["time_start"], event["freq_start"])
        line2 = (event['time_stop'] - event['time_start'],
                 event['freq_stop'] - event['freq_start'])
        line3 = (event['time_stop'], event["freq_stop"])

        lines = [line1, line2, line3]

        orth_slope = -1 / event['slope']

        separation = event["axial_freq"]
        positive_intersections = 0
        negative_intersections = 0

        for i, other_event in field.iterrows():
            # iteration += 1
            positive_intersections = 0
            negative_intersections = 0
            # print('Iteration:', iteration)
            for x, y in lines:
                values = abline(x, y, orth_slope, separation)
                for z in range(2):
                    if is_between(
                        (other_event["time_start"], other_event["freq_start"]),
                        (values[0][z], values[1][z]),
                        (other_event["time_stop"], other_event["freq_stop"])):

                        positive_intersections += 1

                    if is_between(
                        (other_event["time_start"], other_event["freq_start"]),
                        (values[0][z], values[1][z]),
                        (other_event["time_stop"], other_event["freq_stop"])):

                        negative_intersections += 1

        intersections = positive_intersections + negative_intersections

        # print(intersections)

        mainband = -1

        if intersections == 0:
            df.at[index, "Mainband"] = "Mainband"
            # print("Mainband")

        elif intersections == 2:
            if positive_intersections == 1 and negative_intersections == 1:
                df.at[index, "Mainband"] = "Mainband"
                # print("Mainband")
            else:
                df.at[index, "Mainband"] = "Sideband"
                print("Sideband")
        else:
            df.at[index, "Mainband"] = "Unclassified"
            print("Unclassified")

        


toc = time.perf_counter()
print("Time elapsed", toc - tic)

print("Mainbands:", len(df[df["Mainband"] == "Mainband"]))
print("Sidebands:", len(df[df["Mainband"] == "Sideband"]))
print("Unclassified:", len(df[df["Mainband"] == "Unclassified"]))

df.to_csv("dmtracks_bands.csv")
