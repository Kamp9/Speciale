#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "bfpdynamic.cpp"

// using namespace std;

template <typename T>
BFPDynamic<T> gen_bfp(boost::random::mt19937 &rng, const size_t N) {
    vector<T> elems;
    //    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min()+1, numeric_limits<T>::max());
    for(size_t i = 0; i < N; i++){
        elems.push_back(rand_elem(rng));
    }
    boost::random::uniform_int_distribution<T> rand_exp(-10, 10); //rand_exp(numeric_limits<T>::min(), numeric_limits<T>::max());
    BFPDynamic<T> A(elems, rand_exp(rng)); //rand_exp(rng)
    // BFPStatic<T,N> A(elems, 0);
    return A; //BFPStatic<T,N>(A.to_float());
}

template <typename T>
BFPDynamic<T> gen_bfp_pos(boost::random::mt19937 &rng, const size_t N) {
    vector<T> elems;
    //    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
    boost::random::uniform_int_distribution<T> rand_elem(1, numeric_limits<T>::max());
    for(size_t i = 0; i < N; i++){
        elems.push_back(rand_elem(rng));
    }
    boost::random::uniform_int_distribution<T> rand_exp(-10, 10); //rand_exp(numeric_limits<T>::min(), numeric_limits<T>::max());
    BFPDynamic<T> A(elems, rand_exp(rng)); //rand_exp(rng)
    // BFPStatic<T,N> A(elems, 0);
    return A; //BFPStatic<T,N>(A.to_float());
}

template <typename T>
BFPDynamic<T> gen_bfp_neg(boost::random::mt19937 &rng, const size_t N) {
    vector<T> elems;
    //    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min()+1, -1);
    for(size_t i = 0; i < N; i++){
        elems.push_back(rand_elem(rng));
    }
    boost::random::uniform_int_distribution<T> rand_exp(0, 10); //rand_exp(numeric_limits<T>::min(), numeric_limits<T>::max());
    BFPDynamic<T> A(elems, rand_exp(rng)); //rand_exp(rng)
    // BFPStatic<T,N> A(elems, 0);
    return A; //BFPStatic<T,N>(A.to_float());
}

// std::vector<T>(A)

int main(int argc, char *argv[]){
    boost::random::mt19937 rng;
    struct timeval tv;
    gettimeofday(&tv, 0);

    rng.seed(tv.tv_usec);

    auto A = gen_bfp_neg<int8_t>(rng, 5);
    auto B = gen_bfp_neg<int8_t>(rng, 5);
    // auto a = vector<int8_t>{12,108,22,121,73,125,75,84,77,26};
    // auto b = vector<int8_t>{73,122,39,25,98,64,126,73,122,10};

    // auto a = std::vector<int8_t> {-100};
    // auto b = std::vector<int8_t> {-101};
    // BFPDynamic<int8_t> A{a, 6};
    // BFPDynamic<int8_t> B{b, 6};

    // auto a = vector<int8_t>{12,108,22,121,73,125,75,84,77,26};
    // auto b = vector<int8_t>{73,122,39,25,98,64,126,73,122,10};

    // BFPDynamic<int8_t> A{a, -5};
    // BFPDynamic<int8_t> B{b, -8};


    // check_add(A, B);
    check_add(A, B);
    // BFPStatic<int8_t, 1> A100{{-65}, 1};
    // BFPStatic<int8_t, 1> B100{{-65}, 0};



    // cout << bitlog(123) << endl;
    // auto B = gen_bfp_pos<int8_t>(rng, 10);
    // auto B = gen_bfp<int8_t>(rng, 10);
    // vector<int8_t>A2(-15,-107,66,-111,-13);
    // BFPDynamic<int8_t> A{A2, 0, false};
    // cout << bitlog(125) << endl;

    // check_log(A);
    // cout << (A * B).to_float() << endl;
    // cout << A << endl;
    // check_log(A);   
    // cout << floor_log2(int16_t(128)) << endl;
    //     sqrt {{7,63,62,41,109,38,47,43,113,118},9} = 
    // check_sqrt(A);
    // {6,4,4,4,4,4,4,4,4,4}
    // Result of BFP-square root:
    // {{16,80,80,64,112,64,64,64,112,112},1}
    // {{30,90,89,72,118,70,78,74,120,123},1} wanted.

    // cout << A << endl;
    // cout << A.normalized << endl;
    // cout << A.normalize().lazy_list << endl;
    // cout << A.normalized +  << endl;

    return 0;
}
