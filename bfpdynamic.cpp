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

#include <xmmintrin.h>

#include "bfplib.hh"

using namespace std;

template <typename T> ostream &operator<<(ostream &s, const vector<T> &xs){
    s << "{"; for(int i=0;i<xs.size();i++) s << xs[i] << (i+1<xs.size()?"," : ""); s << "}";
    return s;
}

template <typename T, size_t N> ostream &operator<<(ostream &s, const array<T, N> &xs){
    s << "{"; for(int i=0;i<xs.size();i++) s << int(xs[i]) << (i+1<xs.size()?"," : ""); s << "}";
    return s;
}


template <int N> ostream &operator<<(ostream &s, const bitset<N> &bits){
  for(int i=0;i<bits.size();i++) s << bits[i];
  return s;
}

// TODO: Derive double-length type automatically or pass as template parameter.
// typedef int64_t Tx2;
// typedef uint64_t uTx2;

template<typename T> struct uT {typedef T type;};
template<> struct uT<int8_t>  {typedef uint8_t type;};
template<> struct uT<int16_t> {typedef uint16_t type;};
template<> struct uT<int32_t> {typedef uint32_t type;};

template<typename T> struct Tx2 {typedef T type;};
template<> struct Tx2<int8_t>   {typedef int16_t type;};
template<> struct Tx2<int16_t>  {typedef int32_t type;};
template<> struct Tx2<int32_t>  {typedef int64_t type;};

template<typename T> struct uTx2 {typedef T type;};
template<> struct uTx2<int8_t>  {typedef uint16_t type;};
template<> struct uTx2<int16_t> {typedef uint32_t type;};
template<> struct uTx2<int32_t> {typedef uint64_t type;};


// BFPDynamic definition
template <typename T>
struct BFPDynamic: public std::vector<T>{
    int exponent;
    // bool normalized = false;
    // bool lazy = false;
    // std::array<size_t, numeric_limits<T>::digits> lazy_list = {{0}};

    BFPDynamic(int exponent=0) : exponent(exponent) {}
    BFPDynamic(const std::vector<T> &A, int exponent) : std::vector<T>(A), exponent(exponent) {}
    BFPDynamic(const std::vector<double> &V) {
        size_t N = V.size();
        std::vector<T> &A(*this);

        exponent = std::numeric_limits<T>::min();
        for(int i = 0; i < N ;i++){
            int e_V  = floor(log2(fabs(V[i]))) + 1; // TODO: Tjek.
            exponent = std::max(exponent,e_V);
        }
        exponent -= std::numeric_limits<T>::digits;
        double power = pow(2.0,-exponent);
        for(int i=0;i<N;i++){
            assert(round(V[i]*power) >= std::numeric_limits<T>::min());
            assert(round(V[i]*power) <= std::numeric_limits<T>::max());
            // double intpart;
            // double fractpart = modf(V[i]*power, &intpart);
            // Would really like to take care of this in the operators instead
            // But on the other hand, we would need to make calculations to avoid it
            // cout << V[i]*power << endl;
            // if(fractpart == -0.5){
            //     cout << "rounding is happening!" << endl;
            //     A[i] = round(V[i]*power) + 1;
            // }
            // else
            // cout << round(V[i]*power) << endl;
            // cout << V[i]*power << endl;
            A.push_back(round(V[i]*power));
        }
    }

    // BFPDynamic<T> normalize() {
    //     std::vector<T> &A(*this);
    //     const size_t N = this->size();
    //     int min_shifts = numeric_limits<T>::digits + 1;
    //     bool lazy = false;
    //     for(int i=0; i<N; i++){
    //         min_shifts = min(calc_shifts(A[i]), min_shifts);
    //         A[i] <<= min_shifts;
    //         lazy_list[min_shifts] = i;
    //         lazy |= min_shifts != 0;
    //     }
    //     cout << min_shifts << endl;
    //     this->exponent -= min_shifts;
    //     this->normalized = true;
    //     this->lazy = lazy;
    //     return *this;
    // }

    std::vector<double> to_float() const {
        const size_t N = this->size();
        std::vector<double> values(N);
        const std::vector<T> &A(*this);

        for(int i=0;i<N;i++) values[i] = pow(2.0,exponent)*A[i];
        return values;
    }

    friend ostream& operator<<(ostream &s, const BFPDynamic &A) {
        s << "{" << vector<int64_t>(A.begin(),A.end()) << "," << A.exponent << "}";
        return s;
    }
};

// template <typename T,size_t N>
// BFPDynamic<T> operator+(const BFPDynamic<T> &A, const BFPDynamic<T> &B){
//     // Make sure that e_A >= e_B
//     if(A.exponent < B.exponent) return B+A;

//     BFPDynamic<T> AB;
//     bitset<N> carryAB;
//     bitset<N> sAB;
//     int exp_diff = A.exponent - B.exponent;
//     bool carry = false;

//     // for(size_t i=0;i<N;i++){
//     //     Tx2 ABi = A[i] + (B[i] >> exp_diff);
//     //     T   abi = ABi;
//     //     carry  |= (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
//     // }
//     Tx2 max_value = 0;
//     for (size_t i = 0; i < N; i++){
//         Tx2 ABi = Tx2(A[i]) + (B[i] >> exp_diff);
//         max_value = max(max_value, std::abs(ABi) - signbit(ABi));
//     }

//     int shifts = 1 + floor_log2(max_value) - numeric_limits<T>::digits;

//     // if statement on shifts?
//     if(exp_diff == 0){
//         for (size_t i = 0; i < N; i++){
//             AB[i] = Tx2(A[i]) + B[i];
//         }
//         AB.exponent = A.exponent + shifts;
//         return AB;
//     }

//     if(shifts >= 0){
//         for (size_t i = 0; i < N; i++){
//             Tx2 ABi = Tx2(A[i]) + (B[i] >> exp_diff);
//             AB[i] = ABi >> shifts;
//         }
//         AB.exponent = A.exponent + shifts;
//         return AB;
//     }

//     else{
//         shifts = abs(shifts);
//         for (size_t i = 0; i < N; i++){
//             Tx2 ABi = Tx2(A[i]) + (B[i] >> exp_diff);
//             bool rounding = (B[i] >> (exp_diff - 1)) & (!signbit(B[i]));
//             AB[i] = (ABi << shifts) - (rounding << (shifts - 1));
//         }
//         AB.exponent = A.exponent - shifts;
//         return AB;
//     }
// }



template <typename T>
BFPDynamic<T> operator+(const BFPDynamic<T> &A, const BFPDynamic<T> &B){
    if(A.exponent < B.exponent) return B+A;
    size_t N = A.size();

    BFPDynamic<T> AB;

    int exp_diff = A.exponent - B.exponent;
    bool carry = false;

    int exp_diff2 = min(exp_diff, numeric_limits<T>::digits + 1);
    for(size_t i=0;i<N;i++){
        typename Tx2<T>::type ABi = (typename Tx2<T>::type(A[i]) << exp_diff2) + (B[i]);
        ABi = (ABi >> (exp_diff2)) + ((ABi >> (exp_diff2 - 1)) & 1);
        T abi  = ABi;
        carry |= (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
    }
    // Everytime we shift something negative odd number, we have to add one in order to simulate the positive numbers
    if(carry){
        for(size_t i=0;i<N;i++){
            bool sign = signbit((typename Tx2<T>::type(A[i]) << exp_diff2) + B[i]);
            typename Tx2<T>::type ABi = abs((typename Tx2<T>::type(A[i]) << exp_diff2) + B[i]);
            ABi = -sign ^ (ABi >> exp_diff2);
            bool rounding = (ABi & 1);
            AB.push_back((ABi >> 1) + rounding);
        }
    }else{
        for(size_t i=0;i<N;i++){
            typename Tx2<T>::type ABi = (typename Tx2<T>::type(A[i]) << exp_diff2) + (B[i]);
            bool rounding = ((ABi >> (exp_diff2 - 1)) & 1) && (exp_diff2 > 0);
            bool rounding2 = ((ABi & ((1 << exp_diff2) -1)) == (1 << (exp_diff2 - 1))) && signbit(ABi);
            AB.push_back((ABi >> exp_diff2) + rounding - rounding2);
        }
    }
    AB.exponent = A.exponent + carry;

    return AB;
}


            // typename Tx2<T>::type ABi = (typename Tx2<T>::type(A[i]) << exp_diff2) + B[i];
            // typename Tx2<T>::type ABi2 = (typename Tx2<T>::type(abs(A[i])) << exp_diff2) + abs(B[i]);
            // bool rounding = ((ABi2 >> (exp_diff2 - 1)) & 1) && (exp_diff2 > 0);
            // bool rounding2 = ((ABi2 & ((1 << exp_diff2) -1)) == (1 << (exp_diff2 - 1))) && signbit(ABi2);
            // // cout << ABi << endl;
            // AB.push_back(-1 * (((ABi2) >> (exp_diff2))) >> 1);


    // if(carry){
    //     for(size_t i=0;i<N;i++){
    //         cout << "carry" << endl;
    //         typename Tx2<T>::type ABi = (typename Tx2<T>::type(A[i]) << exp_diff2) + (B[i]);
    //         bool rounding = ((ABi >> exp_diff2) & 1);
    //         AB.push_back((ABi >> (exp_diff2 + 1)) + rounding); // + rounding);
    //     }
    // }else{
    //     for(size_t i=0;i<N;i++){
    //         cout << "no carry" << endl;
    //         typename Tx2<T>::type ABi = (typename Tx2<T>::type(A[i]) << exp_diff2) + (B[i]);
    //         bool rounding = ((ABi >> (exp_diff2 - 1)) & 1) && (exp_diff2 > 0);
    //         bool rounding2 = ((ABi & ((1 << exp_diff2) -1)) == (1 << (exp_diff2 - 1))) && signbit(ABi);
    //         AB.push_back((ABi >> exp_diff2) + rounding - rounding2);
    //     }
    // }
    // AB.exponent = A.exponent + carry;

// Working for positive:
// template <typename T>
// BFPDynamic<T> operator+(const BFPDynamic<T> &A, const BFPDynamic<T> &B){
//     if(A.exponent < B.exponent) return B+A;
//     size_t N = A.size();

//     BFPDynamic<T> AB;
//     int exp_diff = A.exponent - B.exponent;
//     bool carry = false;

//     for(size_t i=0;i<N;i++){
//         typename Tx2<T>::type ABi = typename Tx2<T>::type(A[i]) + (B[i] >> exp_diff) + (((B[i] >> (exp_diff - 1)) & 1) && (exp_diff > 0));
//         T abi  = ABi;
//         carry |= (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
//     }

//     AB.exponent = A.exponent + carry;

//     if(carry){
//         for(size_t i=0;i<N;i++){
//             cout << "carry" << endl;
//             typename Tx2<T>::type ABi = typename Tx2<T>::type(A[i]) + (B[i] >> exp_diff);
//             // bool rounding = (((B[i] >> (exp_diff - 2)) & 1) && (exp_diff > 1));
//             bool rounding2 = ABi & 1;
//             cout << ABi << endl;
//             // cout << rounding << endl;
//             cout << rounding2 << endl;
//             AB.push_back((ABi >> 1) +  rounding2); //(rounding || rounding2));
//         }
//     }else{
//         for(size_t i=0;i<N;i++){
//             cout << "no carry" << endl;
//             typename Tx2<T>::type ABi = typename Tx2<T>::type(A[i]) + (B[i] >> exp_diff);
//             bool rounding = (((B[i] >> (exp_diff - 1)) & 1) && (exp_diff > 0));
//             AB.push_back(ABi + rounding);
//         }
//     }
//     AB.exponent = A.exponent + carry;

//     return AB;
// }


// template <typename T>
// BFPDynamic<T> bfp_heateq(const BFPDynamic<T> &A, size_t ydim, size_t xdim){
//     if(A.exponent < B.exponent) return B+A;
    
//     size_t N = A.size();
//     BFPDynamic<T> iteration;

//     max_value = 0;
//     for(size_t y=1; y<ydim-1; y++){
//         for(size_t x=1; x<xdim-1; x++){
//             max_value = (max_value, A[(i-1)*ydim+j] + A[(i+1)*ydim+j] + A[i*ydim+j] + A[i*ydim+j-1] + A[i*ydim+j+1]);
//         }
//     }


//     for(size_t y=1; y<ydim-1; y++){
//         for(size_t x=1; x<xdim-1; x++){
//             iteration.push_back(A[(i-1)*ydim+j] + A[(i+1)*ydim+j] + A[i*ydim+j] + A[i*ydim+j-1] + grid[i*ydim+j+1]);
//         }
//     }

//     return iteration;
// }


template <typename T>
BFPDynamic<T> operator-(const BFPDynamic<T> &A, const BFPDynamic<T> &B){
  size_t N = A.size();
  // assert(N==B.size());
  BFPDynamic<T> nB(N);

  for (size_t i = 0; i < N; i++){ nB.push_back(-B[i]); }
  nB.exponent = B.exponent;

  return A + nB;
}


// template <typename T>
// BFPDynamic<T> operator*(const BFPDynamic<T> &A, const BFPDynamic<T> &B){
//   size_t N = A.size();
//   // assert(N==B.size());
//   BFPDynamic<T> AB;
//   Tx2 max_value = 0;

//   for (size_t i = 0; i < N; i++) {
//     Tx2 ABi = Tx2(A[i]) * B[i];
//     max_value = max(max_value, std::abs(ABi));
//   }
  
//   int shifts = 1 + floor_log2(max_value) - (numeric_limits<T>::digits);
//   if (shifts < 0)
//     shifts = 0;
  
//   for (size_t i = 0; i < N; i++){
//     Tx2 ABi = Tx2(A[i]) * B[i];
//     bool rounding = (ABi >> (shifts - 1)) & 1;
//     bool rounding2 = ((ABi & ((1 << shifts) -1)) == (1 << (shifts - 1))) & signbit(ABi);
//     AB.push_back((ABi >> shifts) + rounding - rounding2);
//   }
//   AB.exponent = A.exponent + B.exponent + shifts;
//   return AB;
// }


// template <typename T>
// BFPDynamic<T> operator/(const BFPDynamic<T> &A, const BFPDynamic<T> &B){
//     size_t N = A.size();
//     // assert(N==B.size());

//     BFPDynamic<T> AB;
//     Tx2 max_value = 0;

//     for (size_t i = 0; i < N; i++) {
//         Tx2 ABi = (Tx2(A[i]) << (numeric_limits<T>::digits + 1)) / B[i];
//         max_value = max(max_value, std::abs(ABi));
//     }
//     int shifts = 1 + floor_log2(max_value) - numeric_limits<T>::digits;
//     if (shifts < 0)
//         shifts = 0;

//     for (size_t i = 0; i < N; i++){
//         Tx2 ABi = (Tx2(A[i]) << (numeric_limits<T>::digits + 1)) / B[i];
//         bool rounding = ((ABi >> (shifts - 1)) & 1);
//         bool rounding2 = ((ABi & ((1 << shifts) -1)) == (1 << (shifts - 1))) & signbit(ABi);  // && needed ?

//         AB.push_back((ABi >> shifts) + rounding - rounding2);
//     }

//     AB.exponent = A.exponent - B.exponent + shifts - numeric_limits<T>::digits - 1;

//     return AB;
// }

// template <typename T>
// BFPDynamic<T> bfp_pow(const BFPDynamic<T> &A, const BFPDynamic<T> &p) {
//     size_t N = A.size();
//     // assert(N==B.size());

//     BFPDynamic<T> Ap;
//     Tx2 max_value = 0;
//     for (size_t i = 0; i < N; i++){
//         Tx2 Api = exp((p[0] << p.exponent) * log(A[i] << A.exponent));
//         max_value = max(max_value, Api);
//     }
//     int shifts = 1 + floor_log2(max_value) - numeric_limits<T>::digits;

//     for (size_t i = 0; i < N; i++){
//         Tx2 Api = exp((p[0] << p.exponent) * log(A[i] << A.exponent));
//         Ap.push_back(Api >> shifts);
//     }

//     Ap.exponent = A.exponent + shifts - numeric_limits<T>::digits;
//     return Ap;
// }

// //Denne bruger floating point
// //https://stackoverflow.com/questions/19611198/finding-square-root-without-using-sqrt-function

// //...nok bedre at basere loesningen paa denne:
// //https://stackoverflow.com/questions/1100090/looking-for-an-efficient-integer-square-root-algorithm-for-arm-thumb2
// template <typename T>
// BFPDynamic<T> bfp_sqrt(const BFPDynamic<T> &A) {
//     size_t N = A.size();
//     // assert(N==B.size());
    
//     BFPDynamic<T> sqrtA;
//     uTx2 max_value = 0;

//     for (size_t i=0; i < N; i++){
//       uTx2 sqrtAi = uTx2(A[i]) << (sizeof(T)*8 + (A.exponent & 1));
//       uTx2 y = (sqrtAi >> 1);
//       uTx2 z = 0;
//       uTx2 y1;

//       while (y != z) {
//         	z = y;
//         	y1 = (y + sqrtAi / y);
//         	y = (y1 >> 1) + (y1 & 1);
//       }
//       max_value = max(max_value, y - (y1 & 1));
//     }

//     int shifts = 1 + floor_log2(max_value) - numeric_limits<T>::digits;

//     if (shifts < 0){
//       shifts = 0;
//     }
//     for (size_t i=0; i < N; i++){
//         uTx2 sqrtAi = uTx2(A[i]) << (numeric_limits<T>::digits + 1 + (A.exponent & 1));
//         uTx2 y = (sqrtAi >> 1);
//         uTx2 z = 0;
//         uTx2 y1;
//         while (y != z) {
//             z = y;
//             y1 = (y + sqrtAi / y);
//             y = (y1 >> 1) + (y1 & 1);
//         }
//         y -= (y1 & 1);

//         sqrtA.push_back((y >> shifts) + ((y >> (shifts - 1) & 1)));
//         }
//     sqrtA.exponent = ((A.exponent + shifts - numeric_limits<T>::digits) >> 1) - signbit(shifts); // (A.exponent & 1);

//     return sqrtA;
// }


// #define BITSPERLONG 32
// #define TOP2BITS(x) ((x & (3L << (BITSPERLONG-2))) >> (BITSPERLONG-2))

// template <typename T>
// BFPDynamic<T> bfp_sqrt2(const BFPDynamic<T> &A) {
//     size_t N = A.size();
//     BFPDynamic<T> sqrtA;
//     uTx2 max_value = 0;

//     for (size_t i=0; i < N; i++){
//         unsigned long x = A[i] << (A.exponent & 1);
//         unsigned long a = 0L;                   /* accumulator      */
//         unsigned long r = 0L;                   /* remainder        */
//         unsigned long e = 0L;                   /* trial product    */

//         for (int j = 0; j < BITSPERLONG; j++){
//             r = (r << 2) + TOP2BITS(x); x <<= 2;
//             a <<= 1;
//             e = (a << 1) + 1;
//             if (r >= e){
//                   r -= e;
//                   a++;
//             }
//         }
//         max_value = max(max_value, a);
//         // max_value = max(max_value, A[i]);
//     }
//     int shifts = 1 + floor_log2(max_value) - numeric_limits<T>::digits;

//     if (shifts > 0){
//         for (size_t i=0; i < N; i++){
//             unsigned long x = A[i] << (A.exponent & 1); // >> (A.exponent & 1);
//             unsigned long a = 0L;                   /* accumulator      */
//             unsigned long r = 0L;                   /* remainder        */
//             unsigned long e = 0L;                   /* trial product    */

//             for (int j = 0; j < BITSPERLONG; j++){
//                 r = (r << 2) + TOP2BITS(x); x <<= 2;
//                 a <<= 1;
//                 e = (a << 1) + 1;
//                 if (r >= e){
//                       r -= e;
//                       a++;
//                 }
//             }
//             bool rounding = ((a >> (shifts - 1)) & 1);
//             sqrtA.push_back((a >> shifts) + rounding);
//         }
//     }else{
//         for (size_t i=0; i < N; i++){
//             unsigned long x = A[i] << (A.exponent & 1); // >> (A.exponent & 1);
//             unsigned long a = 0L;                   /* accumulator      */
//             unsigned long r = 0L;                   /* remainder        */
//             unsigned long e = 0L;                   /* trial product    */

//             for (int j = 0; j < BITSPERLONG; j++){
//                 r = (r << 2) + TOP2BITS(x); x <<= 2;
//                 a <<= 1;
//                 e = (a << 1) + 1;
//                 if (r >= e){
//                       r -= e;
//                       a++;
//                 }
//             }
//             bool rounding = ((a >> (shifts - 1)) & 1);
//             sqrtA.push_back((a << abs(shifts)) + rounding);
//         }
//     }

//     // sqrtA.exponent = 1 + ((A.exponent >> 1) - ((shifts + numeric_limits<T>::digits) >> 1)) + (shifts >> 1);
//     sqrtA.exponent = ((shifts + numeric_limits<T>::digits) >> 1) - (numeric_limits<T>::digits - shifts);

//     return sqrtA;
// }




// // https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
// template <typename T >
// BFPDynamic<T> bfp_invsqrt(const BFPDynamic<T> &A){
//     size_t N = A.size();
//     // assert(N==B.size());

//     BFPDynamic<T> invsqrtA;
//     auto half = BFPDynamic<T>{{1}, -1};
//     auto x2 = A * half;
//     auto threehalfs = BFPDynamic<T>{{3}, -1};
//     auto y = A;

//     y.exponent++;
//     auto i = (BFPDynamic<T>{{0x5f3759df}, 0}) - y;
//     cout << (i[0] >> 20) << endl;
//     i.exponent = i[0] >> 25;

//     auto johni = (i[0] << 8) >> 30;
//     auto res = johni * (threehalfs[0] - (x2[0] * johni * johni));
//     return res;
// }


// template <typename T >
// BFPDynamic<T> bfp_invsqrtfloat(const BFPDynamic<T> &A) {
//     size_t N = A.size();
//     // assert(N==B.size());

//     vector<double> invsqrtAfloat(A.size());
//     for (size_t i=0; i < N; i++){
//         int32_t j;
//         float x2, y;
//         const float threehalfs = 1.5F;
//         x2 = A[i] * 0.5F;

//         y = A[i];
//         j = *(int32_t *) &y;
//         j = 0x5f3759df - (j >> 1);
//         cout << j << endl;
//         y = *(float *) &j;
//         y = y * (threehalfs - (x2 * y * y));
//         invsqrtAfloat.push_back(y);
//     }
//     return BFPDynamic<T>(invsqrtAfloat);
// }


// template <typename T >
// BFPDynamic<T> bfp_exp(const BFPDynamic<T> &A) {
//     BFPDynamic<T> expA;
//     for (size_t i = 0; i < N; i++){
//         double sum = 1.0 + A[i];
//         double term = x;                 // term for k = 1 is just x
//         for (int k = 2; k < 50; k++){
//             term = term * x / (double)k; // term[k] = term[k-1] * x / k
//             sum = sum + term;
//         }
//         expA[i] = sum;
//     }
//     expA.exponent = 0;
//     return expA;
// }


// Vector operations +, -, *, and /
template <typename T> vector<T> operator+(const vector<T> &A, const vector<T> &B){
    assert(A.size() == B.size());
    vector<T> AB(A.size());

    for(int i=0;i<AB.size();i++)
        AB[i] = A[i] + B[i];

    return AB;
}

template <typename T> vector<T> operator-(const vector<T> &A, const vector<T> &B){
    assert(A.size() == B.size());
    vector<T> AB(A.size());

    for(int i=0;i<AB.size();i++)
        AB[i] = A[i] - B[i];

    return AB;
}

template <typename T> vector<T> operator*(const vector<T> &A, const vector<T> &B){
    assert(A.size() == B.size());
    vector<T> AB(A.size());

    for(int i=0;i<AB.size();i++)
        AB[i] = A[i] * B[i];

    return AB;
}

template <typename T> vector<T> operator/(const vector<T> &A, const vector<T> &B){
    assert(A.size() == B.size());
    vector<T> AB(A.size());

    for(int i=0;i<AB.size();i++)
        AB[i] = A[i] / B[i];

    return AB;
}

template <typename T> vector<T> Vpow(const vector<T> &A, const std::vector<T> &p){
    vector<T> Ap(A.size());

    for(int i=0;i<Ap.size();i++)
        Ap[i] = pow(A[i], p[0]);

    return Ap;
}

template <typename T> vector<T> Vsqrt(const vector<T> &A){
    vector<T> sqrtA(A.size());

    for(int i=0;i<sqrtA.size();i++)
        sqrtA[i] = sqrt(A[i]);

    return sqrtA;
}

template <typename T> vector<T> Vinvsqrt(const vector<T> &A){
    vector<T> invsqrtA(A.size());

    for(int i=0;i<invsqrtA.size();i++){
        invsqrtA[i] = 1.0 / std::sqrt(A[i]);
    }

    return invsqrtA;
}


template <typename T> vector<T> Vexp(const vector<T> &A){
    vector<T> expA(A.size());

    for(int i=0;i<expA.size();i++)
        expA[i] = exp(A[i]);

    return expA;
}


// BFPdynamic checkers for operations +, -, *, and /
template <typename T> void check_add(const T& A, const T& B){
    assert(A.size() == B.size());
  size_t N = A.size();
  cout << "******************** Checking addition of: ********************\n"
       << A << " +\n" << B << " = \n\n"
       << A.to_float() << " +\n"  << B.to_float() << "\n\n";

  auto Afloat = A.to_float(), Bfloat = B.to_float();
  auto AB = A+B;
  auto ABfloat = A.to_float() + B.to_float();

  cout << ABfloat << endl;

  vector<bool> sA(N), sB(N);
  for(int i=0;i<N;i++){ sA[i] = signbit(A[i]); sB[i] = signbit(B[i]); }

  cout << "\nFloating point:\n"
       << Afloat << " +\n" << Bfloat << " =\n" << ABfloat << "\n\n";

  cout << "Floating point in BFP-representation:\n"
       << T(Afloat) << " +\n" << T(Bfloat) << " =\n" << T(ABfloat) << "\n\n";

  cout << "Result of BFP-addition:\n"
       << AB << "\n"
       << T(ABfloat) << " wanted.\n\n";

  cout << "Result of BFP-addition converted to floating point:\n"
       << AB.to_float() << "\n"
       << ABfloat << " exact,\n"
       << T(ABfloat).to_float() << " wanted.\n\n";

  cout << "Error compared to exact:\n"
       <<  (AB.to_float() - ABfloat) << "\n"
       << "Error compared to rounded exact:\n"
       << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";
  
  cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");

  // cout << "signs:\n"
  //      << sA << "\n"
  //      << sB << "\n\n";
}

template <typename T> void check_sub(const T& A, const T& B){
  assert(A.size() == B.size());
  size_t N = A.size();
  cout << "******************** Checking subtraction of: ********************\n"
       << A << " -\n" << B << " = \n\n"
       << A.to_float() << " -\n" << B.to_float() << "\n\n";

  auto Afloat = A.to_float(), Bfloat = B.to_float();
  auto AB = A-B;
  auto ABfloat = A.to_float() - B.to_float();

  vector<bool> sA(N), sB(N);
  for(int i=0;i<N;i++){ sA[i] = signbit(A[i]); sB[i] = signbit(B[i]); }

  cout << "\nFloating point:\n"
       << Afloat << " -\n" << Bfloat << " =\n" << ABfloat << "\n\n";

  cout << "Floating point in BFP-representation:\n"
       << T(Afloat) << " -\n" << T(Bfloat) << " =\n" << T(ABfloat) << "\n\n";

  cout << "Result of BFP-subtraction:\n"
       << AB << "\n"
       << T(ABfloat) << " wanted.\n\n";

  cout << "Result of BFP-subtraction converted to floating point:\n"
       << AB.to_float() << "\n"
       << ABfloat << " exact,\n"
       << T(ABfloat).to_float() << " wanted.\n\n";

  cout << "Error compared to exact:\n"
       <<  (AB.to_float() - ABfloat) << "\n"
       << "Error compared to rounded exact:\n"
       << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";
  
  
  cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");

  // cout << "signs:\n"
  //      << sA << "\n"
  //      << sB << "\n\n";
}

template <typename T> void check_mul(const T& A, const T& B){
  assert(A.size() == B.size());
  size_t N = A.size();
  cout << "******************** Checking multiplication of: ********************\n"
       << A << " *\n" << B << " = \n\n"
       << A.to_float() << " *\n"  << B.to_float() << "\n\n";

  auto Afloat = A.to_float(), Bfloat = B.to_float();
  auto AB = A*B;
  auto ABfloat = A.to_float() * B.to_float();

  vector<bool> sA(N), sB(N);
  for(int i=0;i<N;i++){ sA[i] = signbit(A[i]); sB[i] = signbit(B[i]); }

  cout << "\nFloating point:\n"
       << Afloat << " *\n" << Bfloat << " =\n" << ABfloat << "\n\n";

  cout << "Floating point in BFP-representation:\n"
       << T(Afloat) << " *\n" << T(Bfloat) << " =\n" << T(ABfloat) << "\n\n";

  cout << "Result of BFP-multiplication:\n"
       << AB << "\n"
       << T(ABfloat) << " wanted.\n\n";

  cout << "Result of BFP-multiplication converted to floating point:\n"
       << AB.to_float() << "\n"
       << ABfloat << " exact,\n"
       << T(ABfloat).to_float() << " wanted.\n\n";
    
  cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");
  cout << "Error compared to exact:\n"
       <<  (AB.to_float() - ABfloat) << "\n"
       << "Error compared to rounded exact:\n"
       << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";

  // cout << "signs:\n"
  //      << sA << "\n"
  //      << sB << "\n\n";
}

template <typename T> void check_div(const T& A, const T& B){
  assert(A.size() == B.size());
  size_t N = A.size();
  cout << "******************** Checking division of: ********************\n"
       << A << " /\n" << B << " = \n\n"
       << A.to_float() << " /\n"  << B.to_float() << "\n\n";

  auto Afloat = A.to_float(), Bfloat = B.to_float();
  auto AB = A / B;
  auto ABfloat = A.to_float() / B.to_float();

  vector<bool> sA(N), sB(N);
  for(int i=0;i<N;i++){ sA[i] = signbit(A[i]); sB[i] = signbit(B[i]); }

  cout << "\nFloating point:\n"
       << Afloat << " /\n" << Bfloat << " =\n" << ABfloat << "\n\n";

  cout << "Floating point in BFP-representation:\n"
       << T(Afloat) << " /\n" << T(Bfloat) << " =\n" << T(ABfloat) << "\n\n";

  cout << "Result of BFP-division:\n"
       << AB << "\n"
       << T(ABfloat) << " wanted.\n\n";

  cout << "Result of BFP-division converted to floating point:\n"
       << AB.to_float() << "\n"
       << ABfloat << " exact,\n"
       << T(ABfloat).to_float() << " wanted.\n\n";
    
  cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");
  cout << "Error compared to exact:\n"
       <<  (AB.to_float() - ABfloat) << "\n"
       << "Error compared to rounded exact:\n"
       << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";

  // cout << "signs:\n"
  //      << sA << "\n"
  //      << sB << "\n\n";
}

template <typename T> void check_pow(const T& A, const T& p){
  size_t N = A.size();
  cout << "******************** Checking power function of: ********************\n"
       << A << " pow " << p.to_float()[0] << " = \n\n";

  auto Afloat = A.to_float();
  auto Ap = bfp_pow(A, p);
  auto Apfloat = Vpow(A.to_float(), p.to_float());

  cout << "Result of BFP-power:\n"
       << Ap << "\n"
       << T(Apfloat) << " wanted.\n\n";
}


template <typename T> void check_sqrt(const T& A){
  size_t N = A.size();
  cout << "******************** Checking square root function of: ********************\n"
       << "sqrt " << A << " = \n\n";

  // auto Afloat = A.to_float();
  auto sqrtA = bfp_sqrt(A);
  auto sqrtAfloat = Vsqrt(A.to_float());

  cout << "Result of BFP-square root:\n"
       << sqrtA << "\n"
       << T(sqrtAfloat) << " wanted.\n\n";

  cout << "Is the result correct? " << (sqrtA.to_float() == T(sqrtAfloat).to_float()? "Yes.\n" : "No.\n");
}


template <typename T> void check_invsqrt(const T& A){
  size_t N = A.size();
  cout << "******************** Checking inverse square root function of: ********************\n"
       << "invsqrt " << A << " = \n\n";

  // auto Afloat = A.to_float();
  auto invsqrtbfpA = bfp_invsqrt(A);
  auto invsqrtA = bfp_invsqrtfloat(A);
  auto invsqrtAfloat = Vinvsqrt(A.to_float());

  cout << "Result of BFP-inverse square root:\n"
       << invsqrtbfpA << "\n"
       << invsqrtA << "\n"
       << T(invsqrtAfloat) << " wanted.\n\n";

  cout << "Is the result correct? " << (invsqrtA.to_float() == T(invsqrtAfloat).to_float()? "Yes.\n" : "No.\n");
  auto invsqrtA2 = bfp_invsqrtfloat(A);
}
