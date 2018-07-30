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

rects1 = ax.bar(index, means[4*size:] / means[0:size], bar_width,
	        color=sns.xkcd_rgb["denim blue"],
	        error_kw=error_config,
	        alpha=0.9,
	        label='8 bit BFP')

rects2 = ax.bar(index + bar_width, means[4*size:] /  means[size:2*size], bar_width,
	        color=sns.xkcd_rgb["orange"],
	        error_kw=error_config,
	        alpha=0.9,
	        label='16 bit BFP')
	
rects3 = ax.bar(index + 2*bar_width, means[4*size:] /  means[2*size:3*size], bar_width,
	        color=sns.xkcd_rgb["green"],
	        error_kw=error_config,
	        alpha=0.9,
	        label='32 bit BFP')

rects4 = ax.bar(index + 3*bar_width,  means[4*size:] / means[3*size:4*size], bar_width,
	        color=sns.xkcd_rgb["red"],
	        error_kw=error_config,
	        alpha=0.9,
	        label='64 bit BFP')

rects5 = ax.bar(index + 4*bar_width,  native[0:4] , bar_width,
	        color=sns.xkcd_rgb["grey"],
	        error_kw=error_config,
	        alpha=0.9,
	        label='64 bit Floating Point')


ax.plot([-1, 4.5], [1, 1], "k--", linewidth=1.0)

ax.set_ylim(ymin=0.0, ymax=1.5)
ax.set_xlim(xmin=-0.25,xmax=3.75)

plt.title(title)
ax.set_xlabel('Number of elements')
ax.set_ylabel(r'Speedup relative to 64 bit Floating Point') 	

plt.xticks(index + (bar_width * 2), ["10.000.000", "20.000.000", "30.000.000", "40.000.000"])
ax.legend(loc=1)

plt.tight_layout()

plt.savefig('lattice_time_rel.eps', format='eps', dpi=1000)
plt.show()

json.loads
