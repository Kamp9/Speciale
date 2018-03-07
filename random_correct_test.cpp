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
    //    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
    boost::random::uniform_int_distribution<T> rand_elem(0, numeric_limits<T>::max());
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
    boost::random::uniform_int_distribution<T> rand_elem(0, numeric_limits<T>::max());
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

    // BFPStatic<int8_t, 5> C{{32,-56,126,120,-111},0};
    // BFPStatic<int8_t, 5> D{{-83,71,47,98,-72},0};

    // BFPStatic<int8_t, 5> C{{30638,26684,-12140,27759,16812},-10};
    // BFPStatic<int8_t, 5> D{{20393,-9791,-7414,-20592,30398},9};

    // BFPStatic<int8_t, 5> C{{65,-47,-45,27,-69},-6};
    // BFPStatic<int8_t, 5> D{{-88,38,112,4,84},1};

    BFPStatic<int8_t, 5> C{{-32,-119,126,-94,-124},-7};
    BFPStatic<int8_t, 5> D{{32,51,28,-56,59},-1};

    BFPStatic<int8_t, 5> Cpos{{32,119,126,94,124},-7};
    BFPStatic<int8_t, 5> Dpos{{32,51,28,56,59},-1};

    // BFPStatic<int8_t,10> A{{-81, 18, 119, 27, 82, 74, 81, 1, 108, 85},-1};
    // BFPStatic<int8_t,10> B{{-39, -79, 98, -104, 4, 6, 57, 23, 75, 88},-2};
    // cout << (-40.5 + -9.75) << endl;
    // check_add(A,A);
    // check_add(Afp,Afp);
    //    check_add(A,B);
  
    check_add(Cpos,Dpos);
    // check_sub(A,B);
    //    check_mul(A,B);
    // check_div(A,B);

    return 0;
}
