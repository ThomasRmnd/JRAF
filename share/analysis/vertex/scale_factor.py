import argparse
from datetime import datetime

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, nargs="+", help="Input filepaths")
parser.add_argument("--label", type=str, nargs="+", help="Labels of the input files")
args = parser.parse_args()

if not args.input or not args.label:
    print(f"Error: missing required arguments: input or label")
    exit()

if len(args.input) != len(args.label):
    print(f"Error: input and label must have the same length ({len(args.input)} != {len(args.label)}")
    exit()

ninputs = len(args.input)

colors = {
    "P25C": "#6895b9",
    "P25D": "#9b7dbd",
}

fig, ax = plt.subplots(figsize=(10, 5))

for i in range(ninputs):
    timestamps = []
    factors = []
    with open(args.input[i], 'r') as f:
        next(f, None) 
        
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 2: continue
            
            raw_ts, factor = parts
            try:
                if raw_ts.isdigit() or (raw_ts.replace('.', '', 1).isdigit()):
                    ts_obj = datetime.fromtimestamp(float(raw_ts))
                else:
                    ts_obj = datetime.strptime(raw_ts, '%Y-%m-%d %H:%M:%S')
                
                timestamps.append(ts_obj)
                factors.append(float(factor))
            except ValueError as e:
                print(f"Skipping malformed line: {line.strip()} ({e})")

    plot_color = colors.get(args.label[i], "black")
    ax.plot(timestamps, factors, color=plot_color, label=args.label[i], linewidth=1.5)

ax.set_ylabel("Time Correction Factor", fontsize=12)
ax.legend(loc='lower left', frameon=False)

ax.tick_params(direction='in', which='both', top=True, right=True)
ax.minorticks_on()
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))

ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.xticks(rotation=45)

for spine in ax.spines.values():
    spine.set_linewidth(1.2)

fig.tight_layout()
plt.show()