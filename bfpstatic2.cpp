#include <iostream>
#include <vector>
#include <math.h>
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

#include <xmmintrin.h>

#include "bfplib.hh"

using namespace std;

template <typename T> ostream &operator<<(ostream &s, const vector<T> &xs){
    s << "{"; for(int i=0;i<xs.size();i++) s << xs[i] << (i+1<xs.size()?"," : ""); s << "}";
    return s;
}

template <int N> ostream &operator<<(ostream &s, const bitset<N> &bits){
  for(int i=0;i<bits.size();i++) s << bits[i];
  return s;
}

// TODO: Derive double-length type automatically or pass as template parameter.
typedef int64_t Tx2;
typedef uint64_t uTx2;


// BFPStatic definition
template <typename T, size_t N>
struct BFPStatic: public std::array<T,N>{
    int exponent;

  BFPStatic(int exponent=0) : exponent(exponent) {}
  BFPStatic(const std::array<T,N> &A, int exponent) : std::array<T,N>(A), exponent(exponent) {}
  BFPStatic(const std::vector<double> &V) {
        assert(V.size() == N);
    
        std::array<T,N> &A(*this);

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
            A[i] = round(V[i]*power);
        }
    }

    std::vector<double> to_float() const {
        std::vector<double> values(N);
        const std::array<T,N> &A(*this);

        for(int i=0;i<N;i++) values[i] = pow(2.0,exponent)*A[i];
        return values;
    }

    friend ostream& operator<<(ostream &s, const BFPStatic &A) {
        s << "{" << vector<int64_t>(A.begin(),A.end()) << "," << A.exponent << "}";
        return s;
    }
};

// template <typename T,size_t N>
// BFPStatic<T,N> operator+(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
//     // Make sure that e_A >= e_B
//     if(A.exponent < B.exponent) return B+A;

//     BFPStatic<T,N> AB;
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

template <typename T,size_t N>
BFPStatic<T,N> operator+(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    if(A.exponent < B.exponent) return B+A;

    BFPStatic<T,N> AB;
    bitset<N> carryAB;
    bitset<N> sAB;
    int exp_diff = A.exponent - B.exponent;
    bool carry = false;

    for(size_t i=0;i<N;i++){
        Tx2     ABi = A[i] + (B[i] >> exp_diff);
        T       abi = ABi;
        carryAB[i]  = (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
        sAB[i]      = signbit(ABi);
        carry      |= carryAB[i];
    }

    AB.exponent = A.exponent + carry;
    if(carry){
        for(size_t i=0;i<N;i++){
            Tx2 ABi = Tx2(A[i]) + (B[i] >> exp_diff);
            bool rounding = ABi & 1;
            AB[i] = (ABi >> 1) + rounding;
        }
    }
    else{
        for(size_t i=0;i<N;i++){
            Tx2 ABi = Tx2(A[i]) + (B[i] >> exp_diff);
            // does not work when exp_diff = 0
            bool rounding = (B[i] >> (exp_diff - 1)) & 1;
            AB[i] = ABi + rounding;
        }
    }

    AB.exponent = A.exponent + carry;

    return AB;
}


template <typename T,size_t N>
BFPStatic<T,N> operator-(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    BFPStatic<T,N> nB;
    for (size_t i = 0; i < N; i++){ nB[i] = -B[i]; }
    nB.exponent = B.exponent;

    return A + nB;
}


template <typename T,size_t N>
BFPStatic<T,N> operator*(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    BFPStatic<T,N> AB;
    Tx2 max_value = 0;

    for (size_t i = 0; i < N; i++) {
        Tx2 ABi = Tx2(A[i]) * B[i];
        max_value = max(max_value, std::abs(ABi));
    }

    int shifts = 1 + floor_log2(max_value) - (numeric_limits<T>::digits);
    if (shifts < 0)
        shifts = 0;

    for (size_t i = 0; i < N; i++){
        Tx2 ABi = Tx2(A[i]) * B[i];
        bool rounding = (ABi >> (shifts - 1)) & 1;
        bool rounding2 = ((ABi & ((1 << shifts) -1)) == (1 << (shifts - 1))) & signbit(ABi);
        AB[i] = (ABi >> shifts) + rounding - rounding2;
    }
    AB.exponent = A.exponent + B.exponent + shifts;
    return AB;
}


template <typename T,size_t N>
BFPStatic<T,N> operator/(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    BFPStatic<T,N> AB;
    Tx2 max_value = 0;

    for (size_t i = 0; i < N; i++) {
        Tx2 ABi = (Tx2(A[i]) << (numeric_limits<T>::digits + 1)) / B[i];
        max_value = max(max_value, std::abs(ABi));
    }
    int shifts = 1 + floor_log2(max_value) - numeric_limits<T>::digits;
    if (shifts < 0)
        shifts = 0;

    for (size_t i = 0; i < N; i++){
        Tx2 ABi = (Tx2(A[i]) << (numeric_limits<T>::digits + 1)) / B[i];
        bool rounding = ((ABi >> (shifts - 1)) & 1);
        bool rounding2 = ((ABi & ((1 << shifts) -1)) == (1 << (shifts - 1))) & signbit(ABi);  // && needed ?

        AB[i] = (ABi >> shifts) + rounding - rounding2;
    }

    AB.exponent = A.exponent - B.exponent + shifts - numeric_limits<T>::digits - 1;

    return AB;
}

template <typename T, size_t N>
BFPStatic<T, N> bfp_pow(const BFPStatic<T, N> &A, const BFPStatic<T, N> &p) {
    BFPStatic<T,N> Ap;
    Tx2 max_value = 0;
    for (size_t i = 0; i < N; i++){
        Tx2 Api = exp((p[0] << p.exponent) * log(A[i] << A.exponent));
        max_value = max(max_value, Api);
    }
    int shifts = 1 + floor_log2(max_value) - numeric_limits<T>::digits;

    for (size_t i = 0; i < N; i++){
        Tx2 Api = exp((p[0] << p.exponent) * log(A[i] << A.exponent));
        Ap[i] = Api >> shifts;
    }

    Ap.exponent = A.exponent + shifts - numeric_limits<T>::digits;
    return Ap;
}

//https://stackoverflow.com/questions/19611198/finding-square-root-without-using-sqrt-function

template <typename T, size_t N >
BFPStatic<T, N> bfp_sqrt(const BFPStatic<T, N> &A) {
    BFPStatic<T,N> sqrtA;
    uTx2 max_value = 0;

    for (size_t i=0; i < N; i++){
        uTx2 sqrtAi = uTx2(A[i]) << (numeric_limits<T>::digits + 1 + (A.exponent & 1));
        uTx2 y = (sqrtAi >> 1);
        uTx2 z = 0;
        uTx2 y1;

        while (y != z) {
            z = y;
            y1 = (y + sqrtAi / y);
            y = (y1 >> 1) + (y1 & 1);
        }
        max_value = max(max_value, y - (y1 & 1));
    }

    int shifts = 1 + floor_log2(max_value) - numeric_limits<T>::digits;

    if (shifts < 0){
      shifts = 0;
    }
    for (size_t i=0; i < N; i++){
        uTx2 sqrtAi = uTx2(A[i]) << (numeric_limits<T>::digits + 1 + (A.exponent & 1));
        uTx2 y = (sqrtAi >> 1);
        uTx2 z = 0;
        uTx2 y1;
        while (y != z) {
            z = y;
            y1 = (y + sqrtAi / y);
            y = (y1 >> 1) + (y1 & 1);
        }
        y -= (y1 & 1);

        sqrtA[i] = (y >> shifts) + ((y >> (shifts - 1) & 1));

        }
    sqrtA.exponent = ((A.exponent + shifts - numeric_limits<T>::digits) >> 1) - signbit(shifts); // (A.exponent & 1);

    return sqrtA;
}



/**
 * \brief    Fast Square root algorithm, with rounding
 *
 * This does arithmetic rounding of the result. That is, if the real answer
 * would have a fractional part of 0.5 or greater, the result is rounded up to
 * the next integer.
 *      - SquareRootRounded(2) --> 1
 *      - SquareRootRounded(3) --> 2
 *      - SquareRootRounded(4) --> 2
 *      - SquareRootRounded(6) --> 2
 *      - SquareRootRounded(7) --> 3
 *      - SquareRootRounded(8) --> 3
 *      - SquareRootRounded(9) --> 3
 *
 * \param[in] a_nInput - unsigned integer for which to find the square root
 *
 * \return Integer square root of the input value.
 */
// uint32_t SquareRootRounded(uint32_t a_nInput)
// {
//     uint32_t op  = a_nInput;
//     uint32_t res = 0;
//     uint32_t one = 1uL << 30; // The second-to-top bit is set: use 1u << 14 for uint16_t type; use 1uL<<30 for uint32_t type


//     // "one" starts at the highest power of four <= than the argument.
//     while (one > op)
//     {
//         one >>= 2;
//     }

//     while (one != 0)
//     {
//         if (op >= res + one)
//         {
//             op = op - (res + one);
//             res = res +  2 * one;
//         }
//         res >>= 1;
//         one >>= 2;
//     }

//     /* Do arithmetic rounding to nearest integer */
//     if (op > res)
//     {
//         res++;
//     }

//     return res;
// }

// https://stackoverflow.com/questions/1100090/looking-for-an-efficient-integer-square-root-algorithm-for-arm-thumb2

#define BITSPERLONG 32
#define TOP2BITS(x) ((x & (3L << (BITSPERLONG-2))) >> (BITSPERLONG-2))

struct int_sqrt {
      unsigned sqrt,
               frac;
};

template <typename T, size_t N >
BFPStatic<T, N> bfp_sqrt2(const BFPStatic<T, N> &A) {
    BFPStatic<T,N> sqrtA;
    uTx2 max_value = 0;

    for (size_t i=0; i < N; i++){
        unsigned long x = A[i] << (A.exponent & 1);
        unsigned long a = 0L;                   /* accumulator      */
        unsigned long r = 0L;                   /* remainder        */
        unsigned long e = 0L;                   /* trial product    */

        for (int j = 0; j < BITSPERLONG; j++){
            r = (r << 2) + TOP2BITS(x); x <<= 2;
            a <<= 1;
            e = (a << 1) + 1;
            if (r >= e){
                  r -= e;
                  a++;
            }
        }
        max_value = max(max_value, a);
        // max_value = max(max_value, A[i]);
    }
    int shifts = 1 + floor_log2(max_value) - numeric_limits<T>::digits;
    // if (shifts < 0){
    //   shifts = 0;
    // }
    cout << shifts << endl;
    // int shifts = 1 + floor_log2(max_value) - numeric_limits<T>::digits;
    if (shifts > 0){
        for (size_t i=0; i < N; i++){
            unsigned long x = A[i] << (A.exponent & 1); // >> (A.exponent & 1);
            unsigned long a = 0L;                   /* accumulator      */
            unsigned long r = 0L;                   /* remainder        */
            unsigned long e = 0L;                   /* trial product    */

            for (int j = 0; j < BITSPERLONG; j++){
                r = (r << 2) + TOP2BITS(x); x <<= 2;
                a <<= 1;
                e = (a << 1) + 1;
                if (r >= e){
                      r -= e;
                      a++;
                }
            }
            bool rounding = ((a >> (shifts - 1)) & 1);
            sqrtA[i] = (a >> shifts) + rounding;
        }
    }else{
        for (size_t i=0; i < N; i++){
            unsigned long x = A[i] << (A.exponent & 1); // >> (A.exponent & 1);
            unsigned long a = 0L;                   /* accumulator      */
            unsigned long r = 0L;                   /* remainder        */
            unsigned long e = 0L;                   /* trial product    */

            for (int j = 0; j < BITSPERLONG; j++){
                r = (r << 2) + TOP2BITS(x); x <<= 2;
                a <<= 1;
                e = (a << 1) + 1;
                if (r >= e){
                      r -= e;
                      a++;
                }
            }
            bool rounding = ((a >> (shifts - 1)) & 1);
            sqrtA[i] = (a << abs(shifts)) + rounding;
        }
    }


    cout << numeric_limits<T>::digits << endl;
    // sqrtA.exponent = 1 + ((A.exponent >> 1) - ((shifts + numeric_limits<T>::digits) >> 1)) + (shifts >> 1);
    sqrtA.exponent = ((shifts + numeric_limits<T>::digits) >> 1) - (numeric_limits<T>::digits - shifts);

    return sqrtA;
}



// https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
template <typename T, size_t N >
BFPStatic<T, N> bfp_invsqrt(const BFPStatic<T, N> &A){
    BFPStatic<T,N> invsqrtA;
    
    auto half = BFPStatic<T, N>{{1}, -1};
    auto x2 = A * half;
    auto threehalfs = BFPStatic<T, N>{{3}, -1};
    auto y = A;

    y.exponent++;
    auto i = (BFPStatic<T, N>{{0x5f3759df}, 0}) - y;
    
    // cout << (i[0] >> 20) << endl;
    // i.exponent = i[0] >> 25;

    i = i * (threehalfs - (x2 * i * i));
    return i;
}


template <typename T, size_t N >
BFPStatic<T, N> bfp_invsqrtfloat(const BFPStatic<T, N> &A) {
    vector<double> invsqrtAfloat(A.size());
    for (size_t i=0; i < N; i++){
        int32_t j;
        float x2, y;
        const float threehalfs = 1.5F;
        x2 = A[i] * 0.5F;
        y = A[i];
        j = *(int32_t *) &y;
        j = 0x5f3759df - (j >> 1);
        y = *(float *) &j;

        y = y * (threehalfs - (x2 * y * y));
        y = y * (threehalfs - (x2 * y * y));
        invsqrtAfloat[i] = y;
    }
    return BFPStatic<T,N>(invsqrtAfloat);
}


// template <typename T, size_t N >
// BFPStatic<T, N> bfp_exp(const BFPStatic<T, N> &A) {
//     BFPStatic<T,N> expA;
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


// template <typename T> vector<T> Vinvsqrtugly(const vector<T> &A){
//     vector<T> invsqrtA(A.size());

//     __m128 B=_mm_load1_ps(&A); //, B=_mm_load1_ps(&b), C=_mm_load1_ps(&c);
//     // __m128 Thresh=_mm_load1_ps(&thresh);
//     __m128 john = _mm_rsqrt_ss(&B);
//     // invsqrtA[i] = __m128 _mm_rsqrt_ss(__m128 A[i]);

//     return invsqrtA;
// }



template <typename T> vector<T> Vexp(const vector<T> &A){
    vector<T> expA(A.size());

    for(int i=0;i<expA.size();i++)
        expA[i] = exp(A[i]);

    return expA;
}


// BFPstatic checkers for operations +, -, *, and /
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
  auto sqrtA = bfp_sqrt2(A);
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
