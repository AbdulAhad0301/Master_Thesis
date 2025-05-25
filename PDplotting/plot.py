import matplotlib.pyplot as plt

# Sample data from your table
data = {
    70.0: [(3423030, 46.26), (4278788, 26.96), (4706667, 5.41)],
    80.0: [(3423030, 59.96), (4278788, 46.16), (4706667, 36.18)],
    90.0: [(3423030, 72.49), (4278788, 61.29), (4706667, 53.90)],
    100.0: [(3423030, 84.36), (4278788, 74.77), (4706667, 68.69)]
}

plt.figure(figsize=(10, 6))

for P1, points in data.items():
    flows = [x[0] for x in points]
    pressures = [x[1] for x in points]
    plt.plot(flows, pressures, marker='o', label=f'P1 = {P1} bar')

plt.xlabel("Flow Rate (kg/day)")
plt.ylabel("Outlet Pressure P2 (bar)")
plt.title("Outlet Pressure vs Flow Rate for Different Inlet Pressures")
plt.grid(True)
plt.legend(title="Inlet Pressure")
plt.tight_layout()
plt.show()
