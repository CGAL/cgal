import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("build/data/times.txt", skiprows = 1)

x = data[:, 0]
y1 = data[:, 1]
y2 = data[:, 2]



plt.plot(x, y1, 'go--', label='with convex hull')
plt.plot(x, y2, 'ro--', label='without convex hull')
legend = plt.legend(loc='best')
plt.xlabel('#vertices')
plt.ylabel('seconds')

plt.show()