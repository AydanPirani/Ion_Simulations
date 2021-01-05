import os
import pathlib
from datetime import datetime
from data_analysis import gen_lists
from statistics import stdev
from numpy import linspace

# Sets default path to parent SIMION folder
path = pathlib.Path(__file__).parent.absolute()
os.chdir(path.parent.absolute())
print(path)

# Defines output files
log_output = open("Garcia_Lab/log.txt", "a")

# Statistics that will be measured - total and average time
start_time = datetime.now()
total_iterations = 0

# Parameters for the test
samples = 5     # Total amount of "samples", or test points, that will be taken.

ac_start = 0    # Start value for AC and DC voltage ranges
dc_start = 0

ac_end = 10000      # End value for AC and DC voltage ranges
dc_end = 10000

#Standardizing start and end params
ac_interval = (ac_end - ac_start)/samples
dc_interval = (dc_end - dc_start)/samples

ac_start += 0.5 * ac_interval
dc_start += 0.5 * dc_interval

ac_end -= 0.5 * ac_interval
dc_end -= 0.5 * dc_interval

for ac_voltage in linspace(ac_start, ac_end, samples):
    ac_voltage = round(ac_voltage, 3)
    for dc_voltage in linspace(dc_start, dc_end, samples):
        dc_voltage = round(dc_voltage, 3)
        print("AC: {}V, interval = {}V; DC: {}V, interval = {}V".format(ac_voltage, ac_interval, dc_voltage, dc_interval))
        total_iterations += 1

        # Opens Command Line and runs main.lua in SIMION
        print("Connecting to command line...")
        cmd_socket = os.popen("simion --nogui lua Garcia_Lab/main.lua")
        cmd_output = cmd_socket.read()

        print("Connected! Writing output...")
        param_output = open("Garcia_Lab/PAs/current_params.txt", "w")
        param_output.write("{}\n{}".format(ac_voltage,dc_voltage))
        param_output.close()

        # Writes SIMION log to raw data file
        data_output = open("Garcia_Lab/raw_data.txt", "w")
        data_output.write("voltages: ac={}v, dc={}v\n".format(ac_voltage,dc_voltage))
        data_output.write(cmd_output+"\n")
        data_output.close()

        s_x, s_y, e_x, e_y = gen_lists("Garcia_Lab/raw_data.txt")
        efficiency = len(e_x)/10
        concentration_factor = 0
        log_output.write("{},{},{}\n".format(ac_voltage,dc_voltage,efficiency))
        print("{},{},{}\n".format(ac_voltage,dc_voltage,efficiency))

        print("Done with iteration!\n")

end_time = datetime.now()
log_output.close()

# Prints statistics from trial
print("Total time: {}".format(end_time-start_time))
print("Avg. time per iteration: {}".format((end_time-start_time)/total_iterations))
