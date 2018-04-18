import benchpress as bp
# from benchpress.suite_util import BP_ROOT

scripts = [
    ('BFP',  'BFP', [0, 8, 16, 32], [1000000, 10000000, 100000000]),
]

cmd_list = []
for label, name, int_sizes, num_elems in scripts:
    for size in int_sizes:
    	for num in num_elems:
	        full_label = "%s / %s / %s" % (label, size, num)
	        # op_type 1 = addition
	        op_type = 1
	        bash_cmd = "./dynamic_bench_test " + str(size) + " " + str(op_type) + " " + str(num)
	        cmd_list.append(bp.command(bash_cmd, full_label))

# Finally, we build the Benchpress suite, which is written to `--output`
bp.create_suite(cmd_list)

