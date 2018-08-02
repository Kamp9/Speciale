import benchpress as bp
# from benchpress.suite_util import BP_ROOT

# 1111111111111111111111111111111111111111111111

scripts = [
    ('', 'BFP', [8, 16, 32, 64, 0], [10000000, 20000000, 30000000, 40000000]),
]

op_type = 5

cmd_list = []
for label, name, int_sizes, num_elems in scripts:
    for size in int_sizes:
    	for num in num_elems:
	        full_label = "%s / %s / %s" % (op_type, size, num)
	        bash_cmd = "./dynamic_bench_test " + str(size) + " " + str(op_type) + " " + str(num) + " 0"
	        cmd_list.append(bp.command(bash_cmd, full_label))

# Finally, we build the Benchpress suite, which is written to `--output`
bp.create_suite(cmd_list)

