#include "bfplib.hh"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <typeinfo>
#include <time.h>
#include <sys/time.h>
#include <cxxabi.h>


template <typename T> void log2_check(boost::random::mt19937 rng, int N)
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
  }
  fprintf(stderr,"%s check of %d applications of floor_log2<%s>.\n",success?"Successful":"Errors in",N, type_name);
}

int main()
{
  boost::random::mt19937 rng;
  struct timeval tv;
  gettimeofday(&tv, 0);
  rng.seed(tv.tv_usec);
  
  const int N = 1000;

  log2_check<uint8_t> (rng,N);
  log2_check<uint16_t>(rng,N);
  log2_check<uint32_t>(rng,N);
  log2_check<uint64_t>(rng,N);    
    
  return 0;
}
