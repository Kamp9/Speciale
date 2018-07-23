CXX=g++ -O3 -march=native 
#-std=c++11

# all: deterministic_correct_test dynamic_bench_test random_correct_test
all: dynamic_bench_test random_correct_test

# deterministic_correct_test: deterministic_correct_test.cpp
# 	$(CXX) -o deterministic_correct_test deterministic_correct_test.cpp

dynamic_bench_test: bfpdynamic.cpp dynamic_bench_test.cpp
	$(CXX) -o dynamic_bench_test dynamic_bench_test.cpp

dynamic_correct_test: bfpdynamic.cpp dynamic_correct_test.cpp
	$(CXX) -o dynamic_correct_test dynamic_correct_test.cpp

random_correct_test: random_correct_test.cpp
	$(CXX) -o random_correct_test random_correct_test.cpp

%: %.cpp
	$(CXX) -o $@ $^ 

clean:
	-rm -f deterministic_correct_test dynamic_bench_test random_correct_test *.o
