# Data analysis module. Needs to be called separately for each iteration in main.py

import pandas as pd
from io import StringIO

def gen_lists(input_data, radius=5, c_x=131, c_y=131):
    columns = ["Ion N","Events","TOF","Mass","Charge","X","Y","Z","Vt","Vx","Vy","Vz","KE","KE Error"]
    df = pd.DataFrame(columns=columns)

    read_line = False

    lines = open(input_data, "r").read().split("\n")[:-1]

    ac_voltage = 0
    dc_voltage = 0

    for line in lines:
        l = line.split(",")
        if len(l) == len(columns):
            df.loc[len(df)] = list(map(float,l))
        elif "voltages: " in line:
            voltages = line.replace("voltages:","").split(",")
            ac_voltage = voltages[0][4:-1]
            dc_voltage = voltages[1][4:-1]

        if "Begin Next Fly'm" in line:
            read_line = True

    start_x = []
    start_y = []

    end_x = []
    end_y = []

    temp = {"x":0,"y":0}

    print(df)

    for i in range(len(df)):
        data = df.iloc[i,:]

        x = data["X"]
        y = data["Y"]

        if data["TOF"] == 0:
            temp["x"] = x
            temp["y"] = y

        if data["Z"] == 40 and ((x-c_x)**2 + (y-c_y)**2 <= radius**2):
            start_x.append(temp["x"])
            start_y.append(temp["y"])
            end_x.append(x)
            end_y.append(y)

    return start_x, start_y, end_x, end_y
