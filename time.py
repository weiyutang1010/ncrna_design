import numpy as np
import matplotlib.pyplot as plt

with open('./data/eterna/eterna_n50.txt', 'r') as data:
    ids = data.read().split('\n')
    n = [len(line.split(' ')[1]) for line in ids if len(line) > 0]
    ids = [line.split(' ')[0] for line in ids if len(line) > 0]

    time = []
    for p_id in ids:
        with open(f'./results/sampling_time/{p_id}.txt', 'r') as f:
            lines = f.read().split('\n')
            t_line = lines[-2]
            time.append(float(t_line.split(': ')[1]))

# Calculate best-fit line
coefficients = np.polyfit(n, time, 2)
print(coefficients)
poly_line = np.poly1d(coefficients)
x_values = np.linspace(min(n), max(n), 100)
y_values = poly_line(x_values)

plt.scatter(n, time, marker='x', color='orange', label='Exact Time')
plt.plot(x_values, y_values, color='red', label='Best Fit Line: 2.7x^2')

plt.xlabel("Length")
plt.ylabel("Time Taken (s)")
plt.title("Time taken for gradient descent: step=2500, k=1600, b=200")
plt.savefig(f'./time.png', format="png", bbox_inches="tight")

print(n)
print(time)