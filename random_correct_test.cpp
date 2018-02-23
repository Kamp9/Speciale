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
    return BFPStatic<T,N>(elems, 0);
}



int main(){
    boost::random::mt19937 rng;
    rng.seed(time(NULL));

    auto A = gen_bfp<int8_t, 4>(rng);
    auto B = gen_bfp<int8_t, 4>(rng);
    
    BFPStatic<int8_t, 3> C{{-128, -128, -128}, 0};
    BFPStatic<int8_t, 3> D{{-128, -128, -128}, 0};
    auto CD = C+D;
    // auto Q = BFPStatic<int8_t, 3>(C.to_float());
    // cout << Q << endl;
    // check_add(A,A);
    check_add(A,B);
    
    BFPStatic<int8_t, 4> ERR{{-100,73,10,20}, 1};
    cout << BFPStatic<int8_t, 4>(ERR.to_float()).to_float() << endl;
    // check<plus>(A,B);
    // check_add(Afp,Afp);
    // check_multi(A,B);

    return 0;
}