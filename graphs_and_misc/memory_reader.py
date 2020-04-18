import matplotlib
import matplotlib.pyplot as plt
import numpy as np

with open('R_stuff/memory.txt') as memory:
    lines = memory.readlines()
    usage = [int(line.split(' ')[0]) for line in lines if line[0]!='T'][:83]
    usage = [u - usage[0] for u in usage]
    times = [line.strip() for line in lines if line[0]=='T']
    x = [i for i in range(len(usage))][:83]

fig, ax = plt.subplots()

z = np.polyfit(x, usage, 1)
p = np.poly1d(z)


ax.plot(x, p(x), 'r--')
ax.plot(x, usage)

ax.set_xticks([3.04, 6.2+3.04, 24.81+6.2+3.04, 49.81+24.81+6.2+3.04])


ax.set(ylabel='Kilobytes of RAM in Use', xlabel='Minutes Since Start', title='Memory Usage Estimate')
labels = [l.get_text() for l in ax.get_xticklabels()]
for i, g in zip(range(len(labels)), [2, 4, 16, 32]):
    labels[i] = f'End of {g} genome run'


plt.show()

print(sum((3.04, 6.2, 24.81, 49.81)))