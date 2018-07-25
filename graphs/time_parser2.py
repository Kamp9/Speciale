import re
import json

input_file = "my_benchmark"
output_file = input_file + "_time"
parameters_file = output_file + '_parameter'

input_file = open(input_file + ".json", "r")
output_file = open(output_file + ".res", "w")
parameters_file = open(parameters_file + ".res", "w")

pattern = re.compile(r".*elapsed-time: (?P<time>[0-9]+\.[0-9]+)")

result = []
idx = 0
for line in input_file:
    m = re.search(pattern, line)
    if m:
        print m.group('time')

        time = m.group('time')
        if idx % 5 == 0:
            if idx != 0:
                result.append(temp)
            temp = [float(time)]
        else:
            temp += [float(time)]
        idx += 1


result.append(temp)
result_json = json.dumps(result)
output_file.write(result_json)
output_file.close()

