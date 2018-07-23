import benchpress as bp
# from benchpress.suite_util import BP_ROOT

scripts = [
    ('Addition', 'BFP', [8], [50000000, 100000000, 150000000, 200000000], [0, 1]),
]

cmd_list = []
for label, name, int_sizes, num_elems, comp in scripts:
    for size in int_sizes:
    	for c in comp:
	    	for num in num_elems:
		        full_label = "%s / %s / %s / %s" % (label, size, num, c)
		        op_type = 2
		        bash_cmd = "./dynamic_bench_test " + str(size) + " " + str(op_type) + " " + str(num) + " " + str(c)
		        cmd_list.append(bp.command(bash_cmd, full_label))

# Finally, we build the Benchpress suite, which is written to `--output`
bp.create_suite(cmd_list)

