import os
import pathlib
from datetime import datetime
from data_analysis import gen_lists
from statistics import stdev

# Sets default path to parent SIMION folder
path = pathlib.Path(__file__).parent.absolute()
os.chdir(path.parent.absolute())
print(path)

# Defines output files
log_output = open("Garcia_Lab/log.txt", "a")


# Statistics that will be measured - total and average time
start_time = datetime.now()
total_iterations = 0

for ac_voltage in range(0,1):
    for dc_voltage in range(0,1):
        print("AC: {}V, DC: {}V".format(ac_voltage,dc_voltage))
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
