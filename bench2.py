import benchpress as bp
# from benchpress.suite_util import BP_ROOT

# 22222222222222222222222222222222222222222222222

scripts = [
    # ('Addition', 'BFP', [32], [100000000, 200000000, 300000000, 400000000], [0, 1]),
      ('', '', [32], [10000000, 20000000, 30000000, 40000000], [0, 1])
]
op_type = 5

cmd_list = []
for label, name, int_sizes, num_elems, comp in scripts:
    for size in int_sizes:
    	for c in comp:
	    	for num in num_elems:
		        full_label = "%s / %s / %s / %s" % (op_type, size, num, c)
		        bash_cmd = "./dynamic_bench_test " + str(size) + " " + str(op_type) + " " + str(num) + " " + str(c)
		        cmd_list.append(bp.command(bash_cmd, full_label))

# Finally, we build the Benchpress suite, which is written to `--output`
bp.create_suite(cmd_list)

