import numpy as np
import matplotlib.pyplot as plt
import csv
with open("heateq.csv", "rb") as f:
    reader = csv.reader(f)
    heat_matrix = []
    for i in range(50):
	    heat_matrix += [np.array(reader.next()[:-1]).astype(float)]

print heat_matrix

plt.matshow(heat_matrix, cmap='hot')

plt.colorbar()
plt.show()
