import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# List of efficiency curve files
eff_files = [
    "Eff0_65.csv",
    "Eff0_70.csv",
    "Eff0_72.csv",
    "Eff0_74.csv",
    "Eff0_75.csv"
]

# List of speed curve files
speed_files = [
    "SpeedLine1.csv",
    "SpeedLine2.csv",
    "SpeedLine3.csv",
    "SpeedLine4.csv",
    "SpeedLine5.csv",
    "SpeedLine6.csv",
    "SpeedLine7.csv"
]

# Define corresponding RPM labels for speed curves (1500 rpm to 4500 rpm)
speed_rpm = [1500 + i * 500 for i in range(len(speed_files))]

# Initialize the plot
plt.figure(figsize=(10, 7))

# Plot efficiency curves (dashed lines with circle markers)
for file in eff_files:
    df = pd.read_csv(file)
    df.columns = ['x', 'y']
    # Create label from filename, e.g., "Eff0_65.csv" -> "Efficiency 0.65"
    label = file.replace("Eff", "Efficiency ").replace("_", ".").replace(".csv", "")
    plt.plot(df['x'], df['y'], linestyle='--', marker='o', alpha=0.8, color='orange')
    # Optionally, annotate the efficiency curves if desired:
    # mid_idx = len(df) // 2
    # plt.text(df['x'].iloc[mid_idx] + 0.05, df['y'].iloc[mid_idx] + 0.05, label,
    #          fontsize=9, color='orange')

# Plot speed curves (solid lines with "x" markers) and annotate their rightmost point.
for i, file in enumerate(speed_files):
    df = pd.read_csv(file)
    df.columns = ['x', 'y']
    label = f"{speed_rpm[i]} rpm"
    plt.plot(df['x'], df['y'], linestyle='-', marker='x', alpha=0.8, color='brown')
    # Annotate at the rightmost point (max x)
    idx_max = df['x'].idxmax()
    x_last = df['x'].iloc[idx_max]
    y_last = df['y'].iloc[idx_max]
    plt.text(x_last + 0.05, y_last, label, fontsize=9, color='brown')

# Set axis labels and limits
plt.xlabel("Mass flow rate m.sqrt(T)/p")
plt.ylabel("Presure ratio")
plt.xlim(0, 7)
plt.ylim(1, 4)

plt.title("Compressor Map: Efficiency and Speed Curves")
plt.grid(True)
plt.tight_layout()

# Note: Legend is intentionally removed.
plt.show()
