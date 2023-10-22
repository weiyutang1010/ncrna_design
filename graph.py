import sys
import numpy as np
import sys
import matplotlib.pyplot as plt

prob = 0 # 0: - log p(y|x), 1: p(y|x)
if len(sys.argv) >= 2:
    prob = int(sys.argv[1])

i = 0
neg_log_prob = []
prob = []

is_value = False

for line in sys.stdin:
    if i == 0 and len(line) > 0:
        i = 1
        continue

    if line == 'end\n':
        break

    if is_value:
        neg_log_prob.append(float(line))
        prob.append(np.exp(-1 * float(line)))

    if line == 'start\n':
        is_value = True

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

fig, ax1 = plt.subplots()
color = 'red'

ax1.set_xlabel('step')
ax1.set_ylabel('- log p(y | x)', color=color)
ax1.plot(neg_log_prob, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()
color = 'blue'
ax2.set_ylabel('p(y  | x)', color=color)
ax2.plot(prob, color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.title(f'gradient descent')
plt.show()