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

#include <stdio.h>
#include <stdlib.h>

#include "bfpdynamic.cpp"

// Base:
// 	size 
// 	type

// View:
// 	shape
// 	offset
// 	stride
// 	Base*

// Op:
// 	 +
// 	 - 
// 	 *
// 	 /

using namespace std;

// template <typename T, size_t N> ostream &operator<<(ostream &s, const array<T, N> &xs){
//     s << "{"; for(int i=0;i<xs.size();i++) s << int(xs[i]) << (i+1<xs.size()?"," : ""); s << "}";
//     return s;
// }

// template <typename T, size_t N>
// struct Base: public BFPDynamic<T>{
//     Base(BFPDynamic<T> A, int exponent) : BFPDynamic<T>(A) {}
// };



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

    BFPDynamic<T> AB;
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
                AB.push_back((ABi >> 1) + rounding);
            }
        }
    }else{
        for(size_t i=0;i<A.ydim;i++){
            for (size_t j=0;j<A.xdim;j++){
                typename Tx2<T>::type ABi = (typename Tx2<T>::type(A(i,j)) << exp_diff2) + (B(i,j));
                bool rounding = ((ABi >> (exp_diff2 - 1)) & 1) && (exp_diff2 > 0);
                bool rounding2 = ((ABi & ((1 << exp_diff2) -1)) == (1 << (exp_diff2 - 1))) && signbit(ABi);
                AB.push_back((ABi >> exp_diff2) + rounding - rounding2);
            }
        }
    }
    AB.exponent = A.base->exponent + carry;
    return AB;
}


// template <typename T, size_t N>
// double bfp_sum_abs(const BFPStatic<T, N> &A, const size_t ydim, const size_t xdim){
//     BFPStatic<T, N> sumA;
//     // How much do we need???
//     Tx2 accumulator = 0;

//     for(size_t i=1; i<ydim-1; i++){
//         for(size_t j=1; j<xdim-1; j++){
//             accumulator += abs(A[i*ydim+j]);
//         }
//     }
//     return pow(2.0, A.exponent) * accumulator;
// }



// template <typename T>
// BFPDynamic<T> bfp_mul_scalar(const View<T> &A, const double scalar){
//     BFPDynamic<T> Ab;
//     vector<double> V;
//     V.push_back(scalar);

//     BFPDynamic<T> B = BFPDynamic<T>(V);
//     auto Bi = B[0];
//     // should be able to avoid going through A two times with max value for BFP struct in metadata.
//     typename Tx2<T>::type max_value = 0;
//     for(size_t i=0;i<A.ydim;i++){
//         for (size_t j=0;j<A.xdim;j++){
//             typename Tx2<T>::type ABi = typename Tx2<T>::type(A(i,j)) * Bi;
//             max_value = max(max_value, ABi ^ typename Tx2<T>::type(-1 * signbit(ABi)));
//         }
//     }

//     int shifts = 1 + floor_log2(typename uTx2<T>::type(max_value ^ typename Tx2<T>::type(-1 * signbit(max_value)))) - (numeric_limits<T>::digits);
//     if (shifts < 0)
//         shifts = 0;

//     for(size_t i=0;i<A.ydim;i++){
//         for (size_t j=0;j<A.xdim;j++){
//             typename Tx2<T>::type ABi = typename Tx2<T>::type(A(i,j)) * Bi;
//             bool rounding = (ABi >> (shifts - 1)) & 1;
//             bool rounding2 = ((ABi & ((1 << shifts) -1)) == (1 << (shifts - 1))) & signbit(ABi);
//             Ab.push_back((ABi >> shifts) + rounding - rounding2);
//         }
//     }
//     Ab.exponent = A.exponent + B.exponent + shifts;
//     return Ab;
// }


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


extern "C" //This is important to get the function from the -lblas library you will use when compiling
{
    double daxpy_(int *n, double *a, double *A, int *incA, double *B, int *incB);
//The daxpy fortran function shown above multiplies a first matrix 'A' by a constant 'a'
//and adds the result to a second matrix 'B.' Both matrices are of size 'n.'
}

// void daxpy(int n, double a, double *A, int incA, double *B, int incB);
void daxpy(int n, double a, double *A, int incA, double *B, int incB)
{
    daxpy_(&n, &a, A, &incA, B, &incB); //Once again, note the call notation. Important!
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
    // should be able to avoid going through A two times with max value for BFP struct in metadata.
    // typename Tx2<T>::type max_value = 0;
    // for(size_t i=0;i<V.ydim;i++){
    //     for (size_t j=0;j<V.xdim;j++){
    //         // cout << int(V(i,j)) << endl;
    //         typename Tx2<T>::type ABi = typename Tx2<T>::type(V(i,j)) * Bi;
    //         max_value = max(max_value, typename Tx2<T>::type(ABi ^ typename Tx2<T>::type(-signbit(ABi))));
    //         // max_value = max(max_value, ABi ^ typename Tx2<T>::type(-1 * signbit(ABi)));
    //     }
    // }
    // int shifts = 1 + floor_log2(typename uTx2<T>::type(max_value ^ typename Tx2<T>::type(-1 * signbit(max_value)))) - (numeric_limits<T>::digits);
    // if (shifts < 0)
    //     shifts = 0;

    // HMMM
    int shifts = numeric_limits<T>::digits + 1;

    size_t ydim = V.ydim + 2;
    for(size_t i=1;i<V.ydim+1;i++){
        for (size_t j=1;j<V.xdim+1;j++){
            typename Tx2<T>::type ABi = typename Tx2<T>::type(V(i-1,j-1)) * Bi;
            bool rounding = (ABi >> (shifts - 1)) & 1;
            bool rounding2 = ((ABi & ((1 << shifts) -1)) == (1 << (shifts - 1))) & signbit(ABi);
            A[i*ydim+j] = (ABi >> shifts) + rounding - rounding2;
            // Ab.push_back((ABi >> shifts) + rounding - rounding2);
            // cout << int(A[i*ydim+j]) << endl;

        }
    }
    // A.exponent = A.exponent;// + B.exponent + shifts;
}


int main(){
    // A[0] * a + B[0]

    // double A[2] = {2.0, 3.0};
    // double B[2] = {4.0, 4.0};
    // int n = 2, incA=1, incB=1;
    // double a = 51.0;
    // daxpy(n, a, A, incA, B, incB);
    // cout << B[1] << endl;

	// std::array<int8_t, 4> a = {100,7,8,9};
    // std::vector<int8_t> a({2,4,6,8});
	// BFPDynamic<int8_t> A(a, 0);
	// BFPDynamic<int8_t> *pA = &A;

    size_t ydim = 6;
    size_t xdim = 6;

    auto grid_1d = make_grid(ydim, xdim);
    auto A = BFPDynamic<int8_t>(grid_1d);
    auto &Ap = A;

    // xdim??
    auto center = View<int8_t>(&Ap, ydim-2, xdim-2, 1, 1, ydim, 1);
    auto north  = View<int8_t>(&Ap, ydim-2, xdim-2, 0, 1, ydim, 1);
    auto south  = View<int8_t>(&Ap, ydim-2, xdim-2, 2, 1, ydim, 1);
    auto west   = View<int8_t>(&Ap, ydim-2, xdim-2, 1, 0, ydim, 1);
    auto east   = View<int8_t>(&Ap, ydim-2, xdim-2, 1, 2, ydim, 1);

    int iterations = 0;
    int max_iterations = 4;

    double delta = 1.0;
    double epsilon = 0.005; //1e-10;

    // print_grid(A, ydim, xdim);
    // print_view(center);
    // print_view(north);
    // print_view(south);
    // print_view(west);
    // print_view(east);

    while (iterations < max_iterations && delta > epsilon){

        iterations++;

        // auto temp0 = center + north + south + west + east;
        auto bfp_cn    = center + north;

        auto view_cn   = View<int8_t>(&bfp_cn, ydim-2, xdim-2, 0, 0, ydim-2, 1);

        auto bfp_cne   = view_cn + east;
        auto view_cne  = View<int8_t>(&bfp_cne, ydim-2, xdim-2, 0, 0, ydim-2, 1);
        
        auto bfp_sw    = south + west;
        auto view_sw   = View<int8_t>(&bfp_sw, ydim-2, xdim-2, 0, 0, ydim-2, 1);

        auto bfp_cnesw = view_cne + view_sw;
        auto view_cnesw = View<int8_t>(&bfp_cnesw, ydim-2, xdim-2, 0, 0, ydim-2, 1);

        print_view(center);
        bfp_update<int8_t>(A, view_cnesw);

        // center = View<int8_t>(&A, ydim-2, xdim-2, 1, 1, ydim, 1);
        // north  = View<int8_t>(&A, ydim-2, xdim-2, 0, 1, ydim, 1);
        // south  = View<int8_t>(&A, ydim-2, xdim-2, 2, 1, ydim, 1);
        // west   = View<int8_t>(&A, ydim-2, xdim-2, 1, 0, ydim, 1);
        // east   = View<int8_t>(&A, ydim-2, xdim-2, 1, 2, ydim, 1);

        // print_grid(A, ydim, xdim);
    }

	return 0;
}