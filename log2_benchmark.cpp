#include "bfplib.hh"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <typeinfo>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <chrono>

template <typename T> void log2_benchmark(int N, int M, const vector<uint64_t>& testdata)
{
  int Nbits = sizeof(T)*8;

#define clock chrono::high_resolution_clock
  //#define clock chrono::steady_clock

  // Three methods
  // 1) O(n) integer floor_log2_slow(x), 2) FPU (int)floor(log2(x)), and 3) New, fast integer floor_log2(x)
  int msb;
  cerr << "Calculating "<<N<<" calls to "<<Nbits<<"-bit floor_log2_slow(x)... ";
  auto start1 = clock::now();
  for(int i=0;i<N;i++)
    msb = floor_log2_slow<T>(testdata[i%M]);
  auto  end1  = clock::now();
  cerr << chrono::duration<float,milli>(end1-start1).count() << "ms.\n";
  
  volatile int msb2 = msb;
  cerr << "Calculating "<<N<<" calls to "<<Nbits<<"-bit FPU floor(log2(x))... ";
  auto start2 = clock::now();
  for(int i=0;i<N;i++)
    msb = floor(log2(testdata[i%M]));
  auto  end2  = clock::now();
  cerr << chrono::duration<float,milli>(end2-start2).count() << "ms.\n";

  cerr << "Calculating "<<N<<" calls to "<<Nbits<<"-bit new floor_log2(x)... ";
  auto start3 = clock::now();
  for(int i=0;i<N;i++)
    msb = floor_log2<T>(testdata[i%M]);
  auto  end3  = clock::now();
  cerr << chrono::duration<float,milli>(end3-start3).count() << "ms.\n";    
  
  volatile int msb3 = msb;

  cerr << "\n";
}

int main()
{ 
  const int N = 10000000, M = 100;

  // We don't want to measure the random number generation together with the benchmark:  
  auto start0 = clock::now();
  cerr << "Generating random 64-bit numbers... ";
  boost::random::mt19937 rng;
  struct timeval tv;
  gettimeofday(&tv, 0);
  rng.seed(tv.tv_usec);
  
  boost::random::uniform_int_distribution<uint64_t> rand_elem(numeric_limits<uint64_t>::min()+1, numeric_limits<uint64_t>::max());  
  vector<uint64_t> testdata(M);
  for(int i=0;i<M;i++)
    testdata[i] = rand_elem(rng);
  auto  end0  = clock::now();
  cerr << chrono::duration<float,milli>(end0-start0).count() << "ms.\n";

  //  for(auto x: testdata) printf("%lx\n",x);
  
  log2_benchmark<uint8_t> (N,M,testdata);  
  log2_benchmark<uint16_t>(N,M,testdata);
  log2_benchmark<uint32_t>(N,M,testdata);
  log2_benchmark<uint64_t>(N,M,testdata);    
    
  return 0;
}
