#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <ctime>
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


// op types:
// 1 : +
// 2 : -
// 3 : *
// 4 : /
// 5 : sqrt
// 6 : inv sqrt
template <typename T>
void call_op(const int op, const BFPDynamic<T> &A){
     switch(op){
        case 1 :
            A + A;
            break;
        case 2 :
            A - A;
            break;
        case 3 :
            A * A;
            break;
        case 4 :
            A / A;
            break;
        case 5 :
            bfp_sqrt(A);
            break;
        // case 6 :
        //     bfp_invsqrt(A);
        //     break;
        default :
            cout << "Wrong op type! Must be 1, 2, 3, 4, 5, or 6 in second argument" << endl;
    }
}

void call_op_float(const int op, const vector<double> &A){
     switch(op){
        case 1 :
            A + A;
            break;
        case 2 :
            A - A;
            break;
        case 3 :
            A * A;
            break;
        case 4 :
            A / A;
            break;
        case 5 :
            Vsqrt(A);
            break;
        // case 6 :
        //     bfp_invsqrt(A);
        //     break;
        default :
            cout << "Wrong op type! Must be 1, 2, 3, 4, 5, or 6 in second argument" << endl;
    }
}


int main(int argc, char *argv[]){
    boost::random::mt19937 rng;
    struct timeval tv;
    gettimeofday(&tv, 0);
    rng.seed(tv.tv_usec);
    
    int test_type = atoi(argv[1]);
    int op_type   = atoi(argv[2]);
    size_t N      = atoi(argv[3]);

    clock_t begin;
    clock_t end;
    double elapsed_secs;

    // Hmm
    BFPDynamic<int32_t> A;
    vector<double> Afloat;

    BFPDynamic<int8_t> A8;
    BFPDynamic<int16_t> A16;
    BFPDynamic<int32_t> A32;

    switch(test_type){
        case 0 :
            A = gen_bfp<int32_t>(rng, N);
            Afloat = A.to_float();
            begin = clock();

            call_op_float(op_type, Afloat);

            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            break;

        case 8 :
            A8 = gen_bfp<int8_t>(rng, N);
            begin = clock();

            call_op(op_type, A8);

            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            break;

        case 16 :
            A16 = gen_bfp<int16_t>(rng, N);
            begin = clock();
            
            call_op(op_type, A16);
            
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            break;

        case 32 :
            A32 = gen_bfp<int32_t>(rng, N);
            begin = clock();
            
            call_op(op_type, A32);
            
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            break;
        default :
            cout << "Wrong test type! Must be 8, 16, or 32 in first argument." << endl;
    }
    // cout << elapsed_secs << endl;
    cout << "elapsed-time: " << elapsed_secs << endl;
    return 0;
}