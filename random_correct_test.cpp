#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "bfpstatic2.cpp"

// using namespace std;

template <typename T, size_t N>
BFPStatic<T,N> gen_bfp(boost::random::mt19937 &rng) {
    array<T,N> elems;
    //    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min()+1, numeric_limits<T>::max());
    for(size_t i = 0; i < N; i++){
        elems[i] = rand_elem(rng);
    }
    boost::random::uniform_int_distribution<T> rand_exp(-10, 10); //rand_exp(numeric_limits<T>::min(), numeric_limits<T>::max());
    BFPStatic<T,N> A(elems, rand_exp(rng)); //rand_exp(rng)
    // BFPStatic<T,N> A(elems, 0);
    return A; //BFPStatic<T,N>(A.to_float());
}

template <typename T, size_t N>
BFPStatic<T,N> gen_bfp_no_0(boost::random::mt19937 &rng) {
    array<T,N> elems;
    //    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min()+1, numeric_limits<T>::max());
    for(size_t i = 0; i < N; i++){
        auto re = rand_elem(rng);
        while(re == 0)
            re = rand_elem(rng);
        elems[i] = re;
    }
    boost::random::uniform_int_distribution<T> rand_exp(-10, 10); //rand_exp(numeric_limits<T>::min(), numeric_limits<T>::max());
    BFPStatic<T,N> A(elems, rand_exp(rng)); //rand_exp(rng)
    // BFPStatic<T,N> A(elems, 0);
    return A; //BFPStatic<T,N>(A.to_float());
}



template <typename T, size_t N>
BFPStatic<T,N> gen_bfp_pos(boost::random::mt19937 &rng) {
    array<T,N> elems;
    //    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
    boost::random::uniform_int_distribution<T> rand_elem(1, numeric_limits<T>::max());
    for(size_t i = 0; i < N; i++){
        elems[i] = rand_elem(rng);
    }
    boost::random::uniform_int_distribution<T> rand_exp(0, 10); //rand_exp(numeric_limits<T>::min(), numeric_limits<T>::max());
    BFPStatic<T,N> A(elems, rand_exp(rng)); //rand_exp(rng)
    // BFPStatic<T,N> A(elems, 0);
    return A; //BFPStatic<T,N>(A.to_float());
}

template <typename T, size_t N>
BFPStatic<T,N> gen_bfp_neg(boost::random::mt19937 &rng) {
    array<T,N> elems;
    //    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min()+1, -1);
    for(size_t i = 0; i < N; i++){
        elems[i] = rand_elem(rng);
    }
    boost::random::uniform_int_distribution<T> rand_exp(0, 10); //rand_exp(numeric_limits<T>::min(), numeric_limits<T>::max());
    BFPStatic<T,N> A(elems, 0); //rand_exp(rng)
    // BFPStatic<T,N> A(elems, 0);
    return A; //BFPStatic<T,N>(A.to_float());
}

int main(){
    boost::random::mt19937 rng;
    struct timeval tv;
    gettimeofday(&tv, 0);
    rng.seed(tv.tv_usec);

    auto A = gen_bfp_pos<int8_t, 1>(rng);
    BFPStatic<int8_t, 1> p{{2}, 0};
    // auto B = gen_bfp_no_0<int32_t, 1000>(rng);

    // BFPStatic<int8_t, 100> C{{-99,18,49,-101,-42},-8};
    // BFPStatic<int8_t, 100> D{{-15,-107,66,-111,-13},-10};

    // BFPStatic<int8_t, 5> P{{-67,-89,-98,-42,-62},6};
    // BFPStatic<int8_t, 5> Q{{-22,22,95,46,79},10};

    // check_add(A,B);
    check_pow(A, p);

    return 0;
}
