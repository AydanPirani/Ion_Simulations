simion.workbench_program()

io.input("C:\\Users\\aydan\\SIMION\\Garcia_Lab\\PAs\\current_params.txt")

adjustable ac_voltage = tonumber(io.read("*line"))   -- AC voltage
adjustable dc_voltage = tonumber(io.read("*line"))   -- DC voltage
adjustable omega = 1000 * 3.14      -- angular frequency (radians per microsecond)

function segment.fast_adjust()
 adj_elect06 = ac_voltage * (sin(omega * ion_time_of_flight)) + dc_voltage
 adj_elect07 = ac_voltage * (cos(omega * ion_time_of_flight)) + dc_voltage
 adj_elect08 = ac_voltage * (sin(omega * ion_time_of_flight)) + dc_voltage
 adj_elect09 = ac_voltage * (cos(omega * ion_time_of_flight)) + dc_voltage
 adj_elect10 = ac_voltage * (sin(omega * ion_time_of_flight)) + dc_voltage
 adj_elect11 = ac_voltage * (cos(omega * ion_time_of_flight)) + dc_voltage
end
