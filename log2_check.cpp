#include "bfplib.hh"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <typeinfo>
#include <time.h>
#include <sys/time.h>
#include <cxxabi.h>


template <typename T> void test_msb(boost::random::mt19937 rng, int N)
{
  boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min()+1, numeric_limits<T>::max());

  int status;
  char *type_name = abi::__cxa_demangle(typeid(T).name(), 0, 0, &status);
  
  // Tre metoder:
  // 1) O(n) heltals floor_log2(x), 2) (int)floor(log2(x)) (floating point), og 3) Ny, hurtig msb(x)
  bool success = true;
  for(int i=0;i<N;i++){
    T x = rand_elem(rng);
    int msb1 = floor_log2_slow(x);
    int msb2 = (int)floor(log2(x));
    int msb3 = floor_log2(x);

    if(msb1 != msb2 || msb2 != msb3){
      fprintf(stderr,"Error in  floor_log2(%d): !(%d = %d = %d)\n", (int)x, msb1,msb2,msb3);
      success = false;
    }
    //    else fprintf(stderr,"SMUKT! floor_log2(%d): %d = %d = %d\n", (int)x, msb1,msb2,msb3);
  }
  fprintf(stderr,"%s test of floor_log2<%s>.\n",success?"Successful":"Errors in",type_name);
}

int main()
{
  boost::random::mt19937 rng;
  struct timeval tv;
  gettimeofday(&tv, 0);
  rng.seed(tv.tv_usec);
  
  const int N = 1000;

  test_msb<uint8_t> (rng,N);
  test_msb<uint16_t>(rng,N);
  test_msb<uint32_t>(rng,N);
  test_msb<uint64_t>(rng,N);    
    
  return 0;
}
