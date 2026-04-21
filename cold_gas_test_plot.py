# code to plot the cold_gas_test csv file

import matplotlib.pyplot as plt
import csv

filename = "cold_gas_test.csv"

with open(filename, 'r') as f:
    reader = csv.DictReader(f)
    time_s = []
    pressure_bar = []
    temperature_k = []


    for row in reader:
        time_str = row['time_zeroed']
        minutes, seconds = map(float, time_str.split(':'))
        total_seconds = minutes * 60 + seconds
        time_s.append(total_seconds)
        pressure_bar.append(float(row['pressure(bar)']))
        temperature_k.append(float(row['temp(c)']) + 273.15)


print(f"Read {len(time_s)} data points from {filename}.")

# plot pressure and temperature vs time
fig, ax1 = plt.subplots()
color = 'tab:blue'
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Pressure (bar)', color=color)
ax1.plot(time_s, pressure_bar, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel('Temperature (K)', color=color)
ax2.plot(time_s, temperature_k, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()
plt.title('Cold Gas Test: Pressure and Temperature vs Time')
fig.savefig("cold_gas_test_plot.png", dpi=300)
plt.show()