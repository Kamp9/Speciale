#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <bitset>
#include <limits>
#include <array>
#include <vector>
#include <assert.h>
#include <cinttypes>
#include <stdio.h>
#include <functional>
#include <typeinfo>

#include <fstream>
#include <chrono>

#include <stdio.h>
#include <stdlib.h>

#include "bfpdynamic.cpp"


template <typename T>
struct View {
	size_t ydim;
	size_t xdim;
	size_t offset0;
    size_t offset1;
	size_t stride0;    
    size_t stride1;
    BFPDynamic<T> *base;
    // T* base_pts[N];
    View(BFPDynamic<T> *base, size_t ydim, size_t xdim, size_t offset0, size_t offset1, size_t stride0, size_t stride1):
    	base(base), ydim(ydim), xdim(xdim), offset0(offset0), offset1(offset1), stride0(stride0), stride1(stride1){
         //    for(int i = 0; i < ydim; i+= stride){
         //        for(int j = 0; j < xdim; j+= stride){
    	    // 		base_pts[i*xdim + j] = &((*base)[offset + i*xdim + j]);
    	    // 	}
        	// }
        }

    const T operator()(const size_t i, const size_t j) const{
        return (*base)[(i+offset0)*stride0 +
                       (j+offset1)*stride1];
    }
};


// daxpy(l) - Linux man page


// const T operator()(size_t i, size_t j) const{

// }

// T& operator()(size_t i, size_t j) {

// }

// A(i,j)

// for
//     A(i,j);

template <typename T>
BFPDynamic<T> operator+(const View<T> &A, const View<T> &B){
    if(A.base->exponent < B.base->exponent) return B+A;
    assert(A.ydim == B.ydim);
    assert(A.xdim == B.xdim);
    size_t N = A.base->size();

    BFPDynamic<T> AB(N);
    int exp_diff = A.base->exponent - B.base->exponent;
    bool carry = false;

    int exp_diff2 = min(exp_diff, numeric_limits<T>::digits + 1);

    // can we benefit from breaking??
    for(size_t i=0;i<A.ydim;i++){
        for (size_t j=0;j<A.xdim;j++){
            typename Tx2<T>::type ABi = (typename Tx2<T>::type (A(i,j)) << exp_diff2) + (B(i,j));
            ABi = (ABi >> (exp_diff2)) + ((ABi >> (exp_diff2 - 1)) & 1);
            T abi  = ABi;
            carry |= (signbit(A(i,j)) ^ signbit(abi)) & (signbit(B(i,j)) ^ signbit(abi));
        }
    }
    if(carry){
        for(size_t i=0;i<A.ydim;i++){
            for (size_t j=0;j<A.xdim;j++){
                bool sign = signbit((typename Tx2<T>::type(A(i,j)) << exp_diff2) + B(i,j));
                typename Tx2<T>::type ABi = abs((typename Tx2<T>::type(A(i,j)) << exp_diff2) + B(i,j));
                ABi = -sign ^ (ABi >> exp_diff2);
                bool rounding = (ABi & 1);
                AB[i*A.ydim+j] = (ABi >> 1) + rounding;
            }
        }
    }else{
        for(size_t i=0;i<A.ydim;i++){
            for (size_t j=0;j<A.xdim;j++){
                typename Tx2<T>::type ABi = (typename Tx2<T>::type(A(i,j)) << exp_diff2) + (B(i,j));
                bool rounding = ((ABi >> (exp_diff2 - 1)) & 1) && (exp_diff2 > 0);
                bool rounding2 = ((ABi & ((1 << exp_diff2) -1)) == (1 << (exp_diff2 - 1))) && signbit(ABi);
                AB[i*A.ydim+j] = (ABi >> exp_diff2) + rounding - rounding2;
            }
        }
    }
    AB.exponent = A.base->exponent + carry;
    return AB;
}

vector<double> make_grid(const size_t ydim, const size_t xdim){
    double grid[ydim][xdim];

    for (int i=0; i<ydim; i++) {      // Initialize the grid
        for (int j=0;j<xdim;j++) {
            grid[i][j] = 0;
        }
    }

    for (int i=0; i<ydim; i++) {      // And borders
        grid[i][0]      = -273.15;
        grid[i][xdim-1] = -273.15;
    }

    for (int i=0; i<xdim; i++) {
        grid[0][i]      = -273.15;
        grid[ydim-1][i] = 40.0;
    }

    vector<double> grid_1d;
    for (int i = 0; i < ydim; i++){
        for (int j = 0; j < xdim; j++){
            grid_1d.push_back(grid[i][j]);
        }
    }
    return grid_1d;
}
    

template <typename T>
void print_grid(BFPDynamic<T> &A, size_t ydim, size_t xdim){
    cout << "e: " << A.exponent << endl;    
    for(size_t i=0; i<ydim; i++){
        for(size_t j=0; j<xdim; j++){
            printf("%.7f \t", A[i*ydim+j] * pow(2.0, A.exponent));
        }
        cout << endl;
    }
    cout << endl;
}


template <typename T>
void print_view(View<T> &A){
    cout << "e: " << A.base->exponent << endl;
    for(size_t i=0; i<A.ydim; i++){
        for(size_t j=0; j<A.xdim; j++){
            printf("%.7f \t", A(i, j) * pow(2.0, A.base->exponent));
        }
        cout << endl;
    }
    cout << endl;
}


template <typename T>
void bfp_update(BFPDynamic<T> &A, const View<T> &V){ 
    vector<double> b;

    b.push_back(0.2);
    BFPDynamic<T> B = BFPDynamic<T>(b);
    auto Bi = B[0];

    typename Tx2<T>::type max_value = 0;
    for(size_t i=0;i<V.ydim;i++){
        for (size_t j=0;j<V.xdim;j++){
            // cout << int(V(i,j)) << endl;
            typename Tx2<T>::type ABi = typename Tx2<T>::type(V(i,j)) * Bi;
            max_value = max(max_value, typename Tx2<T>::type(ABi ^ typename Tx2<T>::type(-signbit(ABi))));
            // max_value = max(max_value, ABi ^ typename Tx2<T>::type(-1 * signbit(ABi)));
        }
    }
    int shifts = 1 + floor_log2(typename uTx2<T>::type(max_value ^ typename Tx2<T>::type(-1 * signbit(max_value)))) - (numeric_limits<T>::digits);
    if (shifts < 0)
        shifts = 0;

    size_t ydim = V.ydim + 2;
    for(size_t i=1;i<V.ydim+1;i++){
        for (size_t j=1;j<V.xdim+1;j++){
            typename Tx2<T>::type ABi = typename Tx2<T>::type(V(i-1,j-1)) * Bi;
            bool rounding = (ABi >> (shifts - 1)) & 1;
            bool rounding2 = ((ABi & ((1 << shifts) -1)) == (1 << (shifts - 1))) & signbit(ABi);
            A[i*ydim+j] = (ABi >> shifts) + rounding - rounding2;
        }
    }
}


void heat_equation(const size_t ydim, const size_t xdim){

    // remember to change the other one
    // auto grid_1d = make_grid(ydim, xdim);
    auto grid_1d = make_grid(ydim, xdim);

    auto A = BFPDynamic<int8_t>(grid_1d);
    // print_grid(A, ydim, xdim);
    // xdim??
    auto center = View<int8_t>(&A, ydim-2, xdim-2, 1, 1, ydim, 1);
    auto north  = View<int8_t>(&A, ydim-2, xdim-2, 0, 1, ydim, 1);
    auto south  = View<int8_t>(&A, ydim-2, xdim-2, 2, 1, ydim, 1);
    auto west   = View<int8_t>(&A, ydim-2, xdim-2, 1, 0, ydim, 1);
    auto east   = View<int8_t>(&A, ydim-2, xdim-2, 1, 2, ydim, 1);
    typedef std::chrono::high_resolution_clock Clock;
    int iterations = 0;
    int max_iterations = 100;

    double delta = 1.0;
    double epsilon = 0.005; //1e-10;

    ofstream myfile;
    myfile.open ("heateq.csv");

    auto t1 = Clock::now();
    auto t2 = Clock::now();
    auto nano = (t2 - t1).count();

    // check_sin(A);
    t1 = Clock::now();

    while (iterations < max_iterations && delta > epsilon){
        iterations++;
        auto bfp_cn    = center + north;
        auto view_cn   = View<int8_t>(&bfp_cn, ydim-2, xdim-2, 0, 0, ydim-2, 1);

        auto bfp_cne   = view_cn + east;
        auto view_cne  = View<int8_t>(&bfp_cne, ydim-2, xdim-2, 0, 0, ydim-2, 1);
        
        auto bfp_sw    = south + west;
        auto view_sw   = View<int8_t>(&bfp_sw, ydim-2, xdim-2, 0, 0, ydim-2, 1);

        auto bfp_cnesw = view_cne + view_sw;

        auto view_cnesw = View<int8_t>(&bfp_cnesw, ydim-2, xdim-2, 0, 0, ydim-2, 1);

        bfp_update<int8_t>(A, view_cnesw);
    }

    t2 = Clock::now();

    nano = (t2 - t1).count();
    cout << nano << endl;
    for(int i = 0; i < ydim; i++){
        for(int j = 0; j < xdim; j++){
            myfile << (A[i*ydim + j] * pow(2.0, A.exponent)) << ",";
        }
        myfile << endl;
    }
    myfile.close();

}


int main(){

    const size_t ydim = 500;
    const size_t xdim = 500;

    heat_equation(ydim, xdim);
	return 0;
}