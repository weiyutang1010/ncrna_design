import sys
import plotext as plt

i = 0
log = []
is_value = False

for line in sys.stdin:
    if i == 0 and len(line) > 0:
        i = 1
        continue

    if line == 'end\n':
        break

    if is_value:
        log.append(float(line))

    if line == 'start\n':
        is_value = True

# plt.plot(log)
# plt.plot_size(60, 20)
# plt.xlabel('step')

# plt.ylabel('- log p(y|x)')
# plt.title(f'gradient descent')

# plt.show()

plt.plot(log[11:])
plt.plot_size(60, 20)
plt.xlabel('step')
plt.xlim(11)

plt.ylabel('- log p(y|x)')
plt.title(f'gradient descent')

plt.show()