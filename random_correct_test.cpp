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
    boost::random::uniform_int_distribution<T> rand_exp(-10,10); //rand_exp(numeric_limits<T>::min(), numeric_limits<T>::max());
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

    // auto C = gen_bfp_no_0<int8_t, 10>(rng);
    // auto D = gen_bfp_no_0<int8_t, 10>(rng);
    auto A = gen_bfp<int8_t, 10>(rng);
    auto B = gen_bfp<int8_t, 10>(rng);
    check_add(A, B);
    // BFPStatic<int8_t, 3> A{{100, 12, 51}, 0};

    // sqrt {{15,27,19,6,62},9}
    
    // BFPStatic<int32_t, 1> L{{13511351},17};
    // BFPStatic<int8_t, 1> M{{-1},0};

    // Fail in addition?
    // {{107,27,-69,99,68,-21,-122,-47,-38,39,-62,123,-29,-27,26,25,38,51,71,-43,92,-58,-115,-40,5,-50,53,72,118,-34,33,-84,5,65,-17,-63,-28,-77,7,-12},-7} +
    // {{127,60,-125,13,0,107,117,-39,21,-61,44,-34,-52,13,-48,-38,-79,-100,118,-62,38,94,12,119,52,49,-84,-97,-40,-64,-25,-109,37,43,-78,62,47,-70,22,11},0}

    // auto B = gen_bfp_no_0<int32_t, 1000>(rng);

    // BFPStatic<int8_t, 1> A{{123}, 10};
    
    // BFPStatic<int8_t, 1> F{{11}, -2}  ;
    // BFPStatic<int8_t, 1> G{{11}, 1};

    // BFPStatic<int32_t, 1> A{{61},14};
    // BFPStatic<int32_t, 1> A{{61},14};
    // BFPStatic<int8_t, 1> C{{2},10};
    // BFPStatic<int8_t, 100> D{{-15,-107,66,-111,-13},-10};

    // BFPStatic<int8_t, 1> P{{100},0};
    // BFPStatic<int8_t, 1> Q{{100},0};

    // cout << A << endl;
    // cout << bfp_mul_scalar(A, 0.2) << endl;

    // check_sqrt(A);
    // auto S = gen_bfp<int8_t, 20>(rng);
    // auto R = gen_bfp<int8_t, 20>(rng);

    // BFPStatic<int8_t, 1> A100{{-65}, 1};
    // BFPStatic<int8_t, 1> B100{{-65}, 0};

    // check_invsqrt(L);
    // check_sqrt(A);

    return 0;
}
