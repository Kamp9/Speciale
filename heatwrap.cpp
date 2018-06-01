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

#include "bfpstatic2.cpp"

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

template <typename T, size_t N> ostream &operator<<(ostream &s, const array<T, N> &xs){
    s << "{"; for(int i=0;i<xs.size();i++) s << int(xs[i]) << (i+1<xs.size()?"," : ""); s << "}";
    return s;
}

// template <typename T, size_t N>
// struct Base: public BFPStatic<T, N>{
//     Base(BFPStatic<T, N> A, int exponent) : BFPStatic<T, N>(A) {}
// };

template <typename T, size_t N>
struct View : public BFPStatic<T, N> {
	size_t ydim;
	size_t xdim;
	size_t offset;
	size_t stride;
    int exponent;
    T* base_pts[N];
    View(BFPStatic<T, N> *base, size_t ydim, size_t xdim, size_t offset, size_t stride) :
    	BFPStatic<T, N>(*base), ydim(ydim), xdim(xdim), offset(offset), stride(stride) {
            exponent = base->exponent;
            for(int i = 0; i < ydim; i+= stride){
                for(int j = 0; j < xdim; j+= stride){
    	    		base_pts[i*xdim + j] = &((*base)[offset + i*xdim + j]);
    	    	}
        	}
        }
};


template <typename T, size_t N>
BFPStatic<T, N> operator+(const View<T, N> &A, const View<T, N> &B){
    if(A.exponent < B.exponent) return B+A;

    BFPStatic<T, N> AB;
    int exp_diff = A.exponent - B.exponent;
    bool carry = false;

    int exp_diff2 = min(exp_diff, numeric_limits<T>::digits + 1);
    for(size_t i=0;i<N;i++){
        Tx2 ABi = (Tx2(*A.base_pts[i]) << exp_diff2) + (*B.base_pts[i]);
        ABi = (ABi >> (exp_diff2)) + ((ABi >> (exp_diff2 - 1)) & 1);
        T abi  = ABi;
        carry |= (signbit(*A.base_pts[i]) ^ signbit(abi)) & (signbit(*B.base_pts[i]) ^ signbit(abi));
    }
    // Everytime we shift something negative odd number, we have to add one in order to simulate the positive numbers
    if(carry){
        for(size_t i=0;i<N;i++){
            bool sign = signbit((Tx2(*A.base_pts[i]) << exp_diff2) + *B.base_pts[i]);
            Tx2 ABi = abs((Tx2(*A.base_pts[i]) << exp_diff2) + *B.base_pts[i]);
            ABi = -sign ^ (ABi >> exp_diff2);
            bool rounding = (ABi & 1);
            AB[i] = (ABi >> 1) + rounding;
        }
    }else{
        for(size_t i=0;i<N;i++){
            Tx2 ABi = (Tx2(*A.base_pts[i]) << exp_diff2) + (*B.base_pts[i]);
            bool rounding = ((ABi >> (exp_diff2 - 1)) & 1) && (exp_diff2 > 0);
            bool rounding2 = ((ABi & ((1 << exp_diff2) -1)) == (1 << (exp_diff2 - 1))) && signbit(ABi);
            AB[i] = (ABi >> exp_diff2) + rounding - rounding2;
        }
    }
    AB.exponent = A.exponent + carry;
    return AB;
}



int main(){
	std::array<int8_t, 4> a = {100,7,8,9};
	BFPStatic<int8_t, 4> A(a, -51);
	BFPStatic<int8_t, 4> *pA = &A;
    const size_t ydim = 100;
    const size_t xdim = 100;

    size_t stride = 1;
    const size_t slice_N = (ydim-2) * (xdim-2);
    View<int8_t, 4> north = View<int8_t, 4>(pA, 2, 2, 1, stride);
    cout << int(*(north.base_pts[0])) << endl;

 //    View<int8_t, slice_N> north = View<int8_t, slice_N>(pA, ydim, xdim, 1, stride);
 //    View<int8_t, slice_N> south = View<int8_t, slice_N>(pA, ydim, xdim, 1 + 2*ydim, stride);
 //    View<int8_t, slice_N> west = View<int8_t, slice_N>(pA, ydim, xdim, ydim, stride);
 //    View<int8_t, slice_N> east = View<int8_t, slice_N>(pA, ydim, xdim, 2 + ydim, stride);
 //    View<int8_t, slice_N> center = View<int8_t, slice_N>(pA, ydim, xdim, 1 + ydim, stride);

    // cout << north + south + west + east + center << endl;
	return 0;
}