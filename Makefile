all: deterministic_correct_test dynamic_bench_test random_correct_test

deterministic_correct_test: deterministic_correct_test.cpp
	g++ -o deterministic_correct_test deterministic_correct_test.cpp

dynamic_bench_test: dynamic_bench_test.cpp
	g++ -o dynamic_bench_test dynamic_bench_test.cpp

random_correct_test: random_correct_test.cpp
	g++ -o random_correct_test random_correct_test.cpp

clean:
	-rm deterministic_correct_test dynamic_bench_test random_correct_test