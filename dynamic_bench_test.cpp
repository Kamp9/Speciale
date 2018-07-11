#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "bfpdynamic.cpp"


template <typename T>
BFPDynamic<T> gen_bfp(boost::random::mt19937 &rng, const size_t N) {
    vector<T> elems(N);
    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min()+1, numeric_limits<T>::max());
    for(size_t i = 0; i < N; i++){
        auto re = rand_elem(rng);
        while(re == 0)
            re = rand_elem(rng);
        elems[i] = re;
    }
    boost::random::uniform_int_distribution<T> rand_exp(-10, 10);
    BFPDynamic<T> A(elems, rand_exp(rng));
    return A;
}

// template <typename T>
// BFPDynamic<T> gen_bfp(size_t seed, const size_t N) {
//     vector<T> elems(N);
//     srandom(seed);

//     int64_t range_min  = numeric_limits<T>::min();
//     int64_t range_length = int64_t(numeric_limits<T>::max()) - int64_t(numeric_limits<T>::min()) + 1;

//     // cerr << "range min = " << numeric_limits<T>::min() << ", range_max = " << numeric_limits<T>::max() << endl;
//     // cerr << "range length = " << range_length << ", range_min = " << range_min << endl;
    
//     for(size_t i = 0; i < N; i++){
//       long int re = 0;
//       while(re == 0) re = random();
//       elems[i] = (re % range_length) + range_min;
//     }
//     int exponent = (random() % 20) - 10; // -10 til 10
//     BFPDynamic<T> A(elems, exponent);
//     return A;
// }


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
        // case 5 :
        //     bfp_sqrt(A);
        //     break;
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
        // case 5 :
        //     Vsqrt(A);
        //     break;
        // case 6 :
        //     bfp_invsqrt(A);
        //     break;
        default :
            cout << "Wrong op type! Must be 1, 2, 3, 4, 5, or 6 in second argument" << endl;
    }
}


int main(int argc, char *argv[]){


    struct timeval tv;
    gettimeofday(&tv, 0);    
    boost::random::mt19937 rng;    
    // rng.seed(tv.tv_usec);

    rng.seed(42);
    //int rng = 42; 		// Eller brug tv.tv_usec, hvis determinisme ikke oenskes
      
    // We need three arguments
    if (argc != 4){
        cout << "3 arguments needed: test_type, op_type, and size!" << endl;
        return 0;
    }else{
        int test_type = atoi(argv[1]);
        int op_type   = atoi(argv[2]);
        size_t N      = atoi(argv[3]);

        clock_t begin;
        clock_t end;
        double elapsed_secs;

        // Hmm
        BFPDynamic<int32_t> A(N);
<<<<<<< HEAD
        vector<double> Afloat(N);
=======
        vector<double> Afloat;
>>>>>>> f1500f22d0a35d30b58a68257fcb7d2d89060e24

        BFPDynamic<int8_t> A8(N);
        BFPDynamic<int16_t> A16(N);
        BFPDynamic<int32_t> A32(N);

        switch(test_type){
            case 0 :
                A = gen_bfp<int32_t>(rng, N);
                Afloat = A.to_float();
                begin = clock();

                call_op_float(op_type, Afloat);

                end = clock();

		cerr << "Block size: " << (sizeof(decltype(Afloat)::value_type) * Afloat.size()) << endl;		
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                break;

            case 8 :
                A8 = gen_bfp<int8_t>(rng, N);

                begin = clock();
            
                // test.setCallback<TYPE>(f);

                call_op<int8_t>(op_type, A8);

                end = clock();

		cerr << "Block size: " << (sizeof(decltype(A8)::value_type) * A8.size()) << endl;
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                break;

            case 16 :
                A16 = gen_bfp<int16_t>(rng, N);
                begin = clock();
                
                call_op<int16_t>(op_type, A16);
                
                end = clock();

		cerr << "Block size: " << (sizeof(decltype(A16)::value_type) * A16.size()) << endl;		
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                break;

            case 32 :
                A32 = gen_bfp<int32_t>(rng, N);
                begin = clock();
                
                call_op<int32_t>(op_type, A32);
                
                end = clock();

		cerr << "Block size: " << (sizeof(decltype(A32)::value_type) * A32.size()) << endl;				
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                break;
            default :
                cout << "Wrong test type! Must be 8, 16, or 32 in first argument." << endl;
        }
        // cout << elapsed_secs << endl;
        cout << "elapsed-time: " << elapsed_secs << endl;
        return 0;
    }
}
