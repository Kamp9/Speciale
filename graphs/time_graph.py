import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe
import seaborn as sns
import csv

outline=mpe.withStroke(linewidth=1, foreground='black')

plt.rcdefaults()
sns.set(style='ticks',palette='hls')

# parameter_file = "gauss_dw_time_parameter.res"
# parameter_file = file(parameter_file, "r")

native_file = "my_benchmark_time.res"
native_file = file(native_file, "r")

# dw_file = "gauss_dw_time.res"
# dw_file = file(dw_file, "r")

# parameters = json.loads(parameter_file.read())

native = json.loads(native_file.read())
# native = [np.mean(values) for values in native]

# dw = json.loads(dw_file.read())
# dw = [np.mean(values) for values in dw]

title = ""
n_groups = len(native)
fig, ax = plt.subplots()
#diff_arr = native/dw
ones = np.ones(n_groups)
index = np.arange(0,200,50)

bar_width = 0.35
error_config = {'ecolor': '0.3'}

plt.grid()

size = 4

native1 = np.mean(native[0:size], axis=1)
native2 = np.mean(native[size:2*size], axis=1)
native3 = np.mean(native[2*size:3*size], axis=1)
native4 = np.mean(native[3*size:4*size], axis=1)
native5 = np.mean(native[4*size:], axis=1)


rects1 = ax.plot(index, native1, 'o-', alpha=0.9, path_effects=[outline], linewidth=3, color=sns.xkcd_rgb["denim blue"], label='8 bit BFP')
rects2 = ax.plot(index, native2, 'o-', alpha=0.9, path_effects=[outline], linewidth=3, color=sns.xkcd_rgb["orange"], label='16 bit BFP')
rects3 = ax.plot(index, native3, 'o-', alpha=0.9, path_effects=[outline], linewidth=3, color=sns.xkcd_rgb["green"], label='32 bit BFP')
rects4 = ax.plot(index, native4, 'o-', alpha=0.9, path_effects=[outline], linewidth=3, color=sns.xkcd_rgb["red"], label='64 bit BFP')
rects5 = ax.plot(index, native5, 'o-', alpha=0.9, path_effects=[outline], linewidth=3, linestyle='dashdot', color=sns.xkcd_rgb["grey"], label='64 bit Floating Point')

# rects2 = ax.plot(index, dw, alpha=0.9, path_effects=[outline], color=sns.xkcd_rgb["pale red"], label='do_while with iterator')

ax.set_ylim(ymin=0)
ax.set_xlim(xmin=50,xmax=150)

plt.title(title)
ax.set_xlabel('Number of elements')
ax.set_ylabel(r'Total elapsed time (in seconds)')

plt.xticks(index + (bar_width ), ["50.000.000", "100.000.000", "150.000.000", "200.000.000"])
ax.legend(loc=0)

plt.tight_layout()

plt.savefig('gauss_time.eps', format='eps', dpi=1000)
plt.show()

json.loads
