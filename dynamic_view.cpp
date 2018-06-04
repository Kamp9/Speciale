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
struct View : public BFPDynamic<T> {
	size_t ydim;
	size_t xdim;
	size_t offset;
	size_t stride;
    int exponent;
    // T* base_pts[N];
    View(BFPDynamic<T> *base, size_t ydim, size_t xdim, size_t offset, size_t stride) :
    	BFPDynamic<T>(*base), ydim(ydim), xdim(xdim), offset(offset), stride(stride) {
            exponent = base->exponent;
         //    for(int i = 0; i < ydim; i+= stride){
         //        for(int j = 0; j < xdim; j+= stride){
    	    // 		base_pts[i*xdim + j] = &((*base)[offset + i*xdim + j]);
    	    // 	}
        	// }
        }

    const T operator()(const size_t i, const size_t j) const{
        return (*this)[i+offset*stride +
                       j+offset*stride];
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
    if(A.exponent < B.exponent) return B+A;

    assert(A.ydim == B.ydim);
    assert(A.xdim == B.xdim);

    BFPDynamic<T> AB;
    int exp_diff = A.exponent - B.exponent;
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
    AB.exponent = A.exponent + carry;
    return AB;
}


template <typename T, size_t N>
BFPDynamic<T> bfp_mul_scalar(const View<T> &A, const double scalar){
    BFPDynamic<T> Ab;
    vector<double> V;
    V.push_back(scalar);

    BFPDynamic<T> B = BFPDynamic<T>(V);
    auto Bi = B[0];
    // should be able to avoid going through A two times with max value for BFP struct in metadata.
    typename Tx2<T>::type max_value = 0;
    for(size_t i=0;i<A.ydim;i++){
        for (size_t j=0;j<A.xdim;j++){
            typename Tx2<T>::type ABi = typename Tx2<T>::type(A(i,j)) * Bi;
            max_value = max(max_value, ABi ^ typename Tx2<T>::type(-1 * signbit(ABi)));
        }
    }

    int shifts = 1 + floor_log2(typename uTx2<T>::type(max_value ^ typename Tx2<T>::type(-1 * signbit(max_value)))) - (numeric_limits<T>::digits);
    if (shifts < 0)
        shifts = 0;

    for(size_t i=0;i<A.ydim;i++){
        for (size_t j=0;j<A.xdim;j++){
            typename Tx2<T>::type ABi = typename Tx2<T>::type(A(i,j)) * Bi;
            bool rounding = (ABi >> (shifts - 1)) & 1;
            bool rounding2 = ((ABi & ((1 << shifts) -1)) == (1 << (shifts - 1))) & signbit(ABi);
            Ab.push_back((ABi >> shifts) + rounding - rounding2);
        }
    }
    Ab.exponent = A.exponent + B.exponent + shifts;
    return Ab;
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

int main(){
    // A[0] * a + B[0]
    
    double A[2] = {2.0, 3.0};
    double B[2] = {4.0, 4.0};
    int n = 2, incA=1, incB=1;
    double a = 51.0;
    daxpy(n, a, A, incA, B, incB);
    cout << B[1] << endl;

	// std::array<int8_t, 4> a = {100,7,8,9};
    // std::vector<int8_t> a({2,4,6,8});
	// BFPDynamic<int8_t> A(a, 0);
	// BFPDynamic<int8_t> *pA = &A;

    const size_t ydim = 10;
    const size_t xdim = 10;

    auto grid_1d = make_grid(ydim, xdim);
    auto bfp_grid = BFPDynamic<int8_t>(grid_1d);

    BFPDynamic<int8_t> *pA = &bfp_grid;

    size_t stride = 1;

    auto north = View<int8_t>(pA, ydim-2, xdim-2, 1, stride);
    auto south = View<int8_t>(pA, ydim-2, xdim-2, 1 + 2*ydim, stride);
    auto west = View<int8_t>(pA, ydim-2, xdim-2, ydim, stride);
    auto east = View<int8_t>(pA, ydim-2, xdim-2, 2 + ydim, stride);
    auto center = View<int8_t>(pA, ydim-2, xdim-2, 1 + ydim, stride);

    cout << north + south + west + east + center << endl;

	return 0;
}