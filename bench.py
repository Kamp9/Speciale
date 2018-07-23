import benchpress as bp
# from benchpress.suite_util import BP_ROOT

scripts = [
    ('Multiplication', 'BFP', [8, 16, 32, 64, 0], [50000000, 100000000, 150000000, 200000000]),
]

cmd_list = []
for label, name, int_sizes, num_elems in scripts:
    for size in int_sizes:
    	for num in num_elems:
	        full_label = "%s / %s / %s" % (label, size, num)
	        # op_type 1 = addition
	        op_type = 4
	        bash_cmd = "./dynamic_bench_test " + str(size) + " " + str(op_type) + " " + str(num) + " 0"
	        cmd_list.append(bp.command(bash_cmd, full_label))

# Finally, we build the Benchpress suite, which is written to `--output`
bp.create_suite(cmd_list)

