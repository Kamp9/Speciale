#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "bfpdynamic.cpp"

int NUM_REP = 10;

template <typename T>
BFPDynamic<T> gen_bfp(boost::random::mt19937 &rng, const size_t N) {
    vector<T> elems(N);
    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min()+1, numeric_limits<T>::max());
    // boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min()+1, numeric_limits<T>::max());
    for(size_t i = 0; i < N; i++){
        auto re = rand_elem(rng);
        while(re == 0)
           re = rand_elem(rng);
        elems[i] = re;
    }
    boost::random::uniform_int_distribution<T> rand_exp(0, 0);
    BFPDynamic<T> A(elems, rand_exp(rng));
    return A;
}


template <typename T>
BFPDynamic<T> gen_bfp2(boost::random::mt19937 &rng, const size_t N) {
    vector<T> elems(N);
    boost::random::uniform_int_distribution<T> rand_elem(-1, 1);
    for(size_t i = 0; i < N-1; i++){
        elems[i] = rand_elem(rng);
    }
    elems[N-1] = numeric_limits<T>::max();
    boost::random::uniform_int_distribution<T> rand_exp(0, 0);
    BFPDynamic<T> A(elems, rand_exp(rng));
    return A;
}


template <typename T>
BFPDynamic<T> gen_bfp_pos(boost::random::mt19937 &rng, const size_t N) {
    vector<T> elems;
    //    boost::random::uniform_int_distribution<T> rand_elem(numeric_limits<T>::min(), numeric_limits<T>::max());
    boost::random::uniform_int_distribution<T> rand_elem(1, numeric_limits<T>::max());
    for(size_t i = 0; i < N; i++){
        elems.push_back(rand_elem(rng));
    }
    boost::random::uniform_int_distribution<T> rand_exp(0, 0); //rand_exp(numeric_limits<T>::min(), numeric_limits<T>::max());
    BFPDynamic<T> A(elems, rand_exp(rng)); //rand_exp(rng)
    // BFPStatic<T,N> A(elems, 0);
    return A; //BFPStatic<T,N>(A.to_float());
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
T call_op(const int op, const BFPDynamic<T> &A, const BFPDynamic<T> &B){

    clock_t begin;
    clock_t end;
    double elapsed_secs;

    BFPDynamic<T> C(0,0);
    size_t i = random() % A.size();

    begin = clock();

    for(int j = 0; j < NUM_REP; j++){
        switch(op){
            case 1 :
                C = A + A;
                break;
            case 2 :
                C = A - B;
                break;
            case 3 :
                C = A * A;
                break;
            case 4 :
                C = A / A;
                break;
            case 5 :
                C = bfp_sqrt(A);
                break;
            case 6 :
                C = bfp_invsqrt(A);
                break;
            case 7 :
                C = bfp_sin(A);
                break;

            default :
                cout << "Wrong op type! Must be 1, 2, 3, 4, 5, or 6 in second argument" << endl;
        }
    }
    end = clock();
    elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC) /  double(NUM_REP);
    cout << "elapsed-time: " << elapsed_secs << endl;
    return C[i];
}

double call_op_float(const int op, const vector<double> &A, const vector<double> &B){

    vector<double> C;
    size_t i = rand() % A.size();

    clock_t begin;
    clock_t end;
    double elapsed_secs;

    begin = clock();

    for(int j = 0; j < NUM_REP; j++){
        switch(op){
            case 1 :
                C = A + B;
                break;
            case 2 :
                C = A - B;
                break;
            case 3 :
                C = A * B;
                break;
            case 4 :
                C = A / B;
                break;
            case 5 :
                C = Vsqrt(A);
                break;
            case 6 :
                C = Vinvsqrt(A);
                break;
            case 7 :
                C = Vsin(A);
                break;
            default :
                cout << "Wrong op type! Must be 1, 2, 3, 4, 5, or 6 in second argument" << endl;
        }
    }

    end = clock();
    elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC) / double(NUM_REP);

    cout << "elapsed-time: " << elapsed_secs << endl;
    return C[i];
}

template <typename T>
T call_op_comp(const int op, const BFPDynamic<T> &A, const BFPDynamic<T> &B){
    BFPDynamic<T> C(0,0);
    size_t i = rand() % A.size();

     switch(op){
        case 1 :
            // B = A;
            // B.exponent = -20;
            C = plus3(A, B);
            return C[i];
        case 2 :
            C = minus2(A, B);
            return C[i];
        // case 3 :
        //     A * A;
        //     break;
        // case 4 :
        //     A / A;
        //     break;

        case 5 :
            // This is the not optimized one!!!
            C = bfp_sqrt2(A);
            return C[i];
        // case 6 :
        //     bfp_invsqrt(A);
        //     break;
        default :
            cout << "Wrong op type! Must be 1, 2, 3, 4, 5, or 6 in second argument" << endl;
            return C[0];
    }
}

template <typename T>
T call_op_micro(const int op, const BFPDynamic<T> &A){
    BFPDynamic<T> C(A.size());
    size_t i = rand() % A.size();
    switch(op){
        case 1 :
            for(size_t j = 0; j < A.size(); j++){
                C[j] = floor_log2(A[j]);
            };
            return C[i];
        default :
            cout << "Wrong op type! micro bench" << endl;
    }
}




int main(int argc, char *argv[]){
    struct timeval tv;
    gettimeofday(&tv, 0); 

    boost::random::mt19937 rng;
    rng.seed(tv.tv_usec);
    srand(clock());
    // rng.seed(42);
    //int rng = 42; 		// Eller brug tv.tv_usec, hvis determinisme ikke oenskes
        typedef std::chrono::high_resolution_clock Clock;

        int test_type = atoi(argv[1]);
        int op_type   = atoi(argv[2]);
        size_t N      = atoi(argv[3]);
        int comp      = atoi(argv[4]);

        auto t1 = Clock::now();
        auto t2 = Clock::now();
        auto nano = (t2 - t1).count();

        clock_t begin;
        clock_t end;
        double elapsed_secs = 0.0;

        // these should get inside the computations;
        // uint_t?
        BFPDynamic<int8_t>  A8(0, 0);
        BFPDynamic<int16_t> A16(0, 0);
        BFPDynamic<int32_t> A32(0, 0);
        BFPDynamic<int64_t> A64(0, 0);

        BFPDynamic<int8_t>  B8(0,0);
        BFPDynamic<int16_t> B16(0,0);
        BFPDynamic<int32_t> B32(0,0);
        BFPDynamic<int64_t> B64(0,0);

        BFPDynamic<uint8_t>  uA8(0, 0);
        BFPDynamic<uint16_t> uA16(0, 0);
        BFPDynamic<uint32_t> uA32(0, 0);
        BFPDynamic<uint64_t> uA64(0, 0);

        BFPDynamic<uint8_t>  uB8(0,0);
        BFPDynamic<uint16_t> uB16(0,0);
        BFPDynamic<uint32_t> uB32(0,0);
        BFPDynamic<uint64_t> uB64(0,0);

        vector<double> A;
        vector<double> B;
        vector<double> C;

        if(comp == 1){
            switch(test_type){
                case 8:
                    A8 = gen_bfp<int8_t>(rng, N);

                    begin = clock();
                    
                    cout << int(call_op_comp(op_type, A8, A8)) << endl;
                    
                    end = clock();

                    // cerr << "Block size: " << (sizeof(decltype(A32)::value_type) * A32.size()) << endl;              
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    cout << "elapsed-time: " << elapsed_secs << endl;
                    break;  

                case 16:
                    A16 = gen_bfp<int16_t>(rng, N);

                    begin = clock();
                    
                    cout << int(call_op_comp(op_type, A16, A16)) << endl;
                    
                    end = clock();

                    // cerr << "Block size: " << (sizeof(decltype(A32)::value_type) * A32.size()) << endl;              
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    cout << "elapsed-time: " << elapsed_secs << endl;
                    break;  

                case 32:
                    A32 = gen_bfp<int32_t>(rng, N);

                    begin = clock();
                    
                    cout << int(call_op_comp<int32_t>(op_type, A32, A32)) << endl;
                    
                    end = clock();

                    // cerr << "Block size: " << (sizeof(decltype(A32)::value_type) * A32.size()) << endl;              
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    cout << "elapsed-time: " << elapsed_secs << endl;
                    break;  

                case 64:
                    A64 = gen_bfp<int64_t>(rng, N);
 
                    begin = clock();
                    
                    cout << int(call_op_comp<int64_t>(op_type, A64, A64)) << endl;
                    
                    end = clock();

                    // cerr << "Block size: " << (sizeof(decltype(A32)::value_type) * A32.size()) << endl;              
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    cout << "elapsed-time: " << elapsed_secs << endl;
                    break;  
                default :
                    cout << "Wrong test type! Must be 8, 16, or 32 in first argument." << endl;
            }

        }else if (comp == 0){
            switch(test_type){
                case 0 :
                    A = gen_bfp_pos<int64_t>(rng, N).to_float();
                    B = gen_bfp_pos<int64_t>(rng, N).to_float();

                    call_op_float(op_type, A, B);

                    break;

                case 8 :
                    A8 = gen_bfp_pos<int8_t>(rng, N);
                    B8 = gen_bfp_pos<int8_t>(rng, N);

                    call_op(op_type, A8, B8);

                    break;

                case 16 :
                    A16 = gen_bfp_pos<int16_t>(rng, N);
                    B16 = gen_bfp_pos<int16_t>(rng, N);

                    call_op(op_type, A16, B16);
                    
                    break;

                case 32 :
                    A32 = gen_bfp_pos<int32_t>(rng, N);
                    B32 = gen_bfp_pos<int32_t>(rng, N);

                    call_op<int32_t>(op_type, A32, B32);
                    
                    break;

                case 64 :
                    A64 = gen_bfp_pos<int64_t>(rng, N);
                    B64 = gen_bfp_pos<int64_t>(rng, N);
                    call_op<int64_t>(op_type, A64, B64);
                    break;

                default :
                    cout << "Wrong test type! Must be 8, 16, or 32 in first argument." << endl;
            }
        } else if (comp == 2) {
            switch(test_type){
                case 8 :
                    uA8 = gen_bfp_pos<uint8_t>(rng, N);
                    begin = clock();

                    cout << int(call_op_micro(op_type, uA8)) << endl;

                    end = clock();
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    cout << "elapsed-time: " << elapsed_secs << endl;
                    break;
                case 16 :
                    uA16 = gen_bfp_pos<uint16_t>(rng, N);
                    begin = clock();

                    cout << int(call_op_micro(op_type, uA16)) << endl;

                    end = clock();
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    cout << "elapsed-time: " << elapsed_secs << endl;
                    break;
                case 32 :
                    uA32 = gen_bfp_pos<uint32_t>(rng, N);
                    begin = clock();

                    cout << int(call_op_micro(op_type, uA32)) << endl;

                    end = clock();
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    cout << "elapsed-time: " << elapsed_secs << endl;
                    break;
                case 64 :
                    uA64 = gen_bfp_pos<uint64_t>(rng, N);
                    begin = clock();

                    cout << int(call_op_micro(op_type, uA64)) << endl;

                    end = clock();
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    cout << "elapsed-time: " << elapsed_secs << endl;
                    break;
            }
        }
        // cout << elapsed_secs << endl;
        return 0;
    }
// }
