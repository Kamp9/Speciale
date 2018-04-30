#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "bfpdynamic_lazy.cpp"

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
    auto A = gen_bfp_pos<int8_t>(rng, 10);
    // auto B = gen_bfp<int8_t>(rng, 10);
    // auto B = gen_bfp<int8_t>(rng, 10);
    // vector<int8_t>A2(-15,-107,66,-111,-13);
    // BFPDynamic<int8_t> A{A2, 0, false};
    // check_sqrt2(A);
    cout << A << endl;
    check_sqrt(A.normalize());
    // cout << A.normalized << endl;
    // cout << A.normalize().lazy_list << endl;
    // cout << A.normalized +  << endl;

    return 0;
}
