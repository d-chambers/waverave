import numpy as np
import matplotlib.pyplot as plt
import csv


processes=np.arange(1,13,1)

weak_scaling = []
strong_scaling=[]

with open('weak_scaling.out') as txt:
     reader = csv.reader(txt, skipinitialspace=True, delimiter=" ")
     next(reader)
     for row in reader:
         weak_scaling.append(float(row[-2]))
with open('strong_scaling.out') as txt:
     reader = csv.reader(txt, skipinitialspace=True, delimiter=" ")
     next(reader)
     for row in reader:
         strong_scaling.append(float(row[-2]))

fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(10, 5))


ax1.plot(processes, weak_scaling, 'o-')
ax1.set_xlabel("Node Count")
ax1.set_ylabel("Time (s)")
ax1.set_title("Strong Scaling")
ax1.grid(True)
ax1.set_ylim(0, 1.1 * max(weak_scaling))

ax2.plot(processes, strong_scaling, 'o-')
ax2.set_xlabel("Node Count")
ax2.set_ylabel("Time (s)")
ax2.set_title("Weak Scaling")
ax2.grid(True)
ax2.set_ylim(0, 1.1 * max(strong_scaling))
    
plt.savefig("Scaling_results.png")
