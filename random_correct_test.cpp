#include <iostream>
#include <cstdlib>
#include <time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "bfpstatic2.cpp"

using namespace std;

template <typename T, size_t N>
BFPStatic<T,N> gen_bfp(boost::random::mt19937 &rng) {
    array<T,N> elems;
    boost::random::uniform_int_distribution<T> rand(numeric_limits<T>::min(), numeric_limits<T>::max());
    for(size_t i = 0; i < N; i++){
        elems[i] = rand(rng);
    }
    // int exponent = rand(rng) % numeric_limits<T>::digits;

    // assert(V[i]*power >= std::numeric_limits<T>::min());
    // assert(V[i]*power <= std::numeric_limits<T>::max());
    return BFPStatic<T,N>(elems, 0);
}


int main(){
    boost::random::mt19937 rng;
    rng.seed(time(NULL));

    BFPStatic<int8_t, 10> A = gen_bfp<int8_t, 10>(rng);
    BFPStatic<int8_t, 10> B = gen_bfp<int8_t, 10>(rng);
    
    check_add(A,B);
    // check_add(Afp,Afp);
    // check_multi(A,B);

    return 0;
}