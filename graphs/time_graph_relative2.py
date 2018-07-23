import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe
import seaborn as sns
import csv

outline=mpe.withStroke(linewidth=1, foreground='black')

plt.rcdefaults()
sns.set(style='ticks',palette='hls')

# parameter_file = "my_benchmark_time.res"
# parameter_file = file(parameter_file, "r")

native_file = "my_benchmark_time.res"
native_file = file(native_file, "r")

# dw_file = "lattice_dw_time.res"
# dw_file = file(dw_file, "r")

# numpy_file = "lattice_numpy_time.res"
# numpy_file = file(numpy_file, "r")

# parameters = json.loads(parameter_file.read())

native = json.loads(native_file.read())
means = np.array([np.mean(values) for values in native])
native = np.ones(len(native))#[np.mean(values) for values in native]

# dw = json.loads(dw_file.read())
# dw_mean = np.array([np.mean(values) for values in dw])
# dw = native_mean / dw_mean

# numpy = json.loads(numpy_file.read())
# numpy_mean = np.array([np.mean(values) for values in numpy])
# numpy = native_mean / numpy_mean


size = 4
title = ""
n_groups = 4

fig, ax = plt.subplots()
#diff_arr = native/dw
ones = np.ones(n_groups)
index = np.arange(0,4,1)

bar_width = 0.15
error_config = {'ecolor': '0.3'}

plt.grid()
# print(numpy_mean)
# print(dw_mean)


# native2 = np.mean(means[size:2*size], axis=1)
# native3 = np.mean(means[2*size:3*size], axis=1)
# native4 = np.mean(means[3*size:4*size], axis=1)
# native5 = np.mean(means[4*size:], axis=1)

# print means[0:size]
# print means[4*size:]
# print means[4*size:] / means[0:size]
	
rects1 = ax.bar(index, native[0:4], bar_width,
	        color=sns.xkcd_rgb["black"],
	        error_kw=error_config,
	        alpha=0.9,
	        label='64 bit BFP')

rects2 = ax.bar(index + bar_width,  means[0:size] /  means[size:2*size], bar_width,
	        color=sns.xkcd_rgb["lightish purple"],
	        error_kw=error_config,
	        alpha=0.9,
	        label='64 bit BFP - optimized')

ax.plot([-1, 4], [1, 1], "k--", linewidth=1.0)

ax.set_ylim(ymin=0.0, ymax=max(means) + 0.5)
ax.set_xlim(xmin=-0.25,xmax=3.5)

plt.title(title)
ax.set_xlabel('Number of elements')
ax.set_ylabel(r'Speedup relative to 64 bit Floating Point') 	

plt.xticks(index + (bar_width*0.5), ["50.000.000", "100.000.000", "150.000.000", "200.000.000"])
ax.legend(loc=1)

plt.tight_layout()

plt.show()

json.loads
