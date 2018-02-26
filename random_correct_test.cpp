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
    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
    for(size_t i = 0; i < N; i++){
        elems[i] = rand_elem(rng);
    }
    boost::random::uniform_int_distribution<T> rand_exp(-10, 10); //rand_exp(numeric_limits<T>::min(), numeric_limits<T>::max());
    BFPStatic<T,N> A(elems, 0); //rand_exp(rng)
    // BFPStatic<T,N> A(elems, 0);
    return A; //BFPStatic<T,N>(A.to_float());
}



int main(){
    boost::random::mt19937 rng;
    rng.seed(time(NULL));

    auto A = gen_bfp<int8_t, 4>(rng);
    auto B = gen_bfp<int8_t, 4>(rng);

    BFPStatic<int8_t, 3> C{{-64, 72, 80}, 0};
    BFPStatic<int8_t, 3> D{{40, -80, 120}, 0};
    // 2400 5600 9600

    // auto CD = C+D;
    // auto Q = BFPStatic<int8_t, 3>(C.to_float());
    // cout << Q << endl;
    // check_add(A,A);
    check_mul(C,D);
    
    // check<plus>(A,B);
    // check_add(Afp,Afp);
    // check_multi(A,B);

    return 0;
}