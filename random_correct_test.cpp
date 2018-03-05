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
    BFPStatic<T,N> A(elems, rand_exp(rng)); //rand_exp(rng)
    // BFPStatic<T,N> A(elems, 0);
    return A; //BFPStatic<T,N>(A.to_float());
}

template <typename T, size_t N>
BFPStatic<T,N> gen_bfp_no_0(boost::random::mt19937 &rng) {
    array<T,N> elems;
    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
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


int main(){
    boost::random::mt19937 rng;
    rng.seed(time(NULL));

    auto A = gen_bfp_no_0<int8_t, 5>(rng);
    auto B = gen_bfp_no_0<int8_t, 5>(rng);

    // BFPStatic<int8_t, 5> C{{65,-47,-45,27,-69},-6};
    // BFPStatic<int8_t, 5> D{{-88,38,112,4,84},1};
    // cout << D.to_float() << endl;
    // BFPStatic<int16_t, 1> E{{-31179},1};
    // BFPStatic<int16_t, 1> F{{30007},8};

    // vector<double> v{-203161600};
    // BFPStatic<int8_t, 1> V = BFPStatic<int8_t, 1>(v);
    // cout << V << endl;
    // auto CD = C+D;
    // auto Q = BFPStatic<int8_t, 3>(C.to_float());
    // cout << Q << endl;
    // check_add(A,A);
    check_div(A,B);

    // + 0 5 5 6 6 6 7 7
    // - 0 1 2 2 3 4 5 5 5
    // check<plus>(A,B);
    // check_add(Afp,Afp);
    // check_multi(A,B);

    return 0;
}