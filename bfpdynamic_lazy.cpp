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
// typedef uint16_t uTx2;

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
    int lazy_min = numeric_limits<T>::digits; // numeric_limits<T>::digits + 2; // 1?
    int normalized = false;
    int contains_nagative = false;
    // std::array<int64_t, numeric_limits<T>::digits> lazy_list;
    std::vector<int64_t> lazy_list;

    BFPDynamic(int exponent=0) : exponent(exponent) {}

    BFPDynamic(std::vector<T> &A, int exponent) : std::vector<T>(A), exponent(exponent) {
        // lazy_list.fill(-1);
    }

    BFPDynamic(std::vector<double> &V) {
        size_t N = V.size();
        std::vector<T> &A(*this);

        exponent = std::numeric_limits<T>::min();
        for(int i = 0; i < N ;i++){
            int e_V  = floor(log2(fabs(V[i]))) + 1;
            exponent = std::max(exponent,e_V);
        }

        exponent -= std::numeric_limits<T>::digits;
        double power = pow(2.0,-exponent);
        for(int i=0;i<N;i++){
            assert(round(V[i]*power) >= std::numeric_limits<T>::min());
            assert(round(V[i]*power) <= std::numeric_limits<T>::max());
            A.push_back(round(V[i]*power));
        }
    }

    BFPDynamic<T> normalize() {
        std::vector<T> &A(*this);
        const size_t N = this->size();
        int min_shifts = numeric_limits<int>::max();
        int bits = numeric_limits<T>::digits + 1;
        for(int i=0; i<N; i++){
            // min_shifts = min(calc_shifts(bits, A[i]), min_shifts);
            typename uT<T>::type absA = A[i] ^ (typename uT<T>::type(-1 * signbit(A[i])));
            min_shifts = min(bits - floor_log2(absA) - 2, min_shifts);
            A[i] <<= min_shifts;
            lazy_list.push_back(min_shifts);
            contains_nagative |= signbit(A[i]);
            // lazy_list[min_shifts] = i;
        }
        exponent -= min_shifts;
        lazy_min = min_shifts;
        normalized = true;
        return *this;
    }

    BFPDynamic<T> denormalize() {
        cout << *this << endl;

        std::vector<T> &A(*this);
        const size_t N = this->size();
        cout << lazy_list << endl << endl;;
        for(int i=0; i<N; i++){
            if (lazy_list[i] - lazy_min > 0) {
                bool rounding = (A[i] >> ((lazy_list[i] - lazy_min) - 1)) & 1;
                // cout << ((A[i] >> ((lazy_list[i] - lazy_min) - 1)) & 1) << endl;
                A[i] = (A[i] >> (lazy_list[i] - lazy_min)) + rounding;
            }else{
                A[i] = (A[i] >> (lazy_list[i] - lazy_min));
            } 
        }
        normalized = false;
        return *this;
    }

    std::vector<double> to_float() {
        if(normalized){
            this->denormalize();
        }

        const size_t N = this->size();
        std::vector<double> values(N);
        const std::vector<T> &A(*this);

        for(int i=0;i<N;i++){
            values[i] = pow(2.0,exponent) * A[i];
        }
        return values;
    }



    friend ostream& operator<<(ostream &s, const BFPDynamic &A) {
        s << "{" << vector<int64_t>(A.begin(),A.end()) << "," << A.exponent << "}";
        return s;
    }
};


template <typename T>
BFPDynamic<T> operator+(BFPDynamic<T> &A, BFPDynamic<T> &B){
    if (!A.normalized){
        A.normalize();
    }

    if (!B.normalized){
        B.normalize();
    }
    // We have to normalize first in case they change the magnitude of A.exponent or B.exponent
    if(A.exponent < B.exponent) return B+A;

    size_t N = A.size();
    BFPDynamic<T> AB;

    int bits = numeric_limits<T>::digits + 1;
    int min_shifts = numeric_limits<int>::max();
    int exp_diff = A.exponent - B.exponent;
    bool carry = false;
    for(size_t i=0;i<N;i++){
        T Ai = A[i] >> (A.lazy_list[i] - A.lazy_min);
        T Bi = B[i] >> (B.lazy_list[i] - B.lazy_min);
        typename Tx2<T>::type ABi = typename Tx2<T>::type(Ai) + (Bi >> exp_diff);

        typename uTx2<T>::type absABi = ABi ^ (typename Tx2<T>::type(-1 * signbit(ABi)));
        min_shifts = min(bits - floor_log2(absABi) - 1, min_shifts);

        bool rounding = false;
        if (exp_diff > 0){
            rounding = (Bi >> (exp_diff - 1)) & 1;
        }

        bool rounding2 = (ABi & 1) & (min_shifts == 0);
        ABi += rounding | rounding2;
        cout << (rounding | rounding2) << endl;
        // cout << rounding << endl;
        // cout << rounding2 << endl << endl;

        AB.lazy_list.push_back(min_shifts);
        AB.push_back((ABi << min_shifts) >> 1);
    }
    AB.exponent = A.exponent - min_shifts + 1;
    AB.lazy_min = min_shifts;
    AB.normalized = true;
    return AB;
}


    // size_t N = A.size();
    // assert(N == B.size());


    // int exp_diff = A.exponent - B.exponent;
    // bool carry = false;
    // for(size_t i=0;i<N;i++){
    //     Tx2 ABi = Tx2(A[i]) + (B[i] >> (exp_diff));
    //     T abi = ABi;
    //     carry |= (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
    //     AB.lazy_list[carry] = i;

    //     AB.push_back(ABi >> carry);
    // }
    // // We should make a case when exp_diff are 0

    // AB.exponent = A.exponent + carry;
    // // AB.lazy = carry;
    // // AB.normalized = true;
    // return AB;


template <typename T>
BFPDynamic<T> operator-(const BFPDynamic<T> &A, const BFPDynamic<T> &B){
  size_t N = A.size();
  // assert(N==B.size());
  BFPDynamic<T> nB(N);

  for (size_t i = 0; i < N; i++){ nB.push_back(-B[i]); }
  nB.exponent = B.exponent;

  return A + nB;
}


template <typename T>
BFPDynamic<T> operator*(BFPDynamic<T> &A, BFPDynamic<T> &B){
    BFPDynamic<T> AB;
    size_t N = A.size();
    int min_shifts = numeric_limits<int>::max();
    int bits = numeric_limits<T>::digits + 1;
    int double_bits = bits * 2;

    if (!A.normalized){
        A.normalize();
    }
    if (!B.normalized){
        B.normalize();
    }

    for (size_t i = 0; i < N; i++){
        typename Tx2<T>::type ABi = typename Tx2<T>::type(A[i] >> (A.lazy_list[i] - A.lazy_min)) * (B[i] >> (B.lazy_list[i] - B.lazy_min));
        min_shifts = min(calc_shifts(double_bits, ABi), min_shifts);
        AB.lazy_list.push_back(min_shifts);

        bool rounding = (ABi & 1) & (min_shifts == 0);
        ABi += rounding;

        AB.push_back((ABi << min_shifts) >> bits);
    }
    AB.lazy_min = min_shifts;
    AB.normalized = true;
    AB.exponent = A.exponent + B.exponent - min_shifts + bits;
    return AB;
}


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

//Denne bruger floating point
//https://stackoverflow.com/questions/19611198/finding-square-root-without-using-sqrt-function

//...nok bedre at basere loesningen paa denne:
//https://stackoverflow.com/questions/1100090/looking-for-an-efficient-integer-square-root-algorithm-for-arm-thumb2
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
//             z = y;
//             y1 = (y + sqrtAi / y);
//             y = (y1 >> 1) + (y1 & 1);
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


template <typename T>
BFPDynamic<T> bfp_sqrt(BFPDynamic<T> &A) {
    BFPDynamic<T> sqrtA;
    size_t N = A.size();

    int min_shifts = numeric_limits<int>::max();

    int bitsdouble = (numeric_limits<T>::digits + 1) * 2;
    int bits = bitsdouble >> 1;
    if (!A.normalized){
        A.normalize();
    }

    typename uTx2<T>::type a;
    typename uTx2<T>::type root;
    typename uTx2<T>::type rem;

    for (size_t i=0; i < N; i++){
        a = typename uTx2<T>::type(A[i] >> (A.lazy_list[i] - A.lazy_min)) << (bits - 2 + (A.exponent & 1));
        rem = 0;
        root = 0;
        for(int i=0; i<bits; i++){
            root <<= 1;
            rem = ((rem << 2) + (a >> (bitsdouble - 2)));
            a <<= 2;
            root++;
            if(root <= rem){
                rem -= root;
                root++;
            }else
                root--;
        }
        root >>= 1;

        bool rounding = root < rem;
        root += rounding;
        min_shifts = min(bits - floor_log2(root) - 2, min_shifts);
        sqrtA.lazy_list.push_back(min_shifts);
        sqrtA.push_back((root << min_shifts));
    }
    sqrtA.exponent = 1 + (A.exponent >> 1) - min_shifts - (bits >> 1);
    sqrtA.lazy_min = min_shifts;
    sqrtA.normalized = true;
    return sqrtA;
}

unsigned short bfp_sqrt3 (unsigned a){
    // a = a >> 2;
    unsigned rem = 0;
    unsigned root = 0;
        for(int i=0; i<16; i++){
            root <<= 1;
            rem = ((rem << 2) + (a >> 30));
            a <<= 2;
            root ++;
            if(root <= rem){
                rem -= root;
                root++;
            }
            else
            root--;
        }
    return (unsigned short)(root >> 1);
}

// https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
template <typename T >
BFPDynamic<T> bfp_invsqrt(const BFPDynamic<T> &A){
    size_t N = A.size();
    // assert(N==B.size());

    BFPDynamic<T> invsqrtA;
    auto half = BFPDynamic<T>{{1}, -1};
    auto x2 = A * half;
    auto threehalfs = BFPDynamic<T>{{3}, -1};
    auto y = A;

    y.exponent++;
    auto i = (BFPDynamic<T>{{0x5f3759df}, 0}) - y;
    cout << (i[0] >> 20) << endl;
    i.exponent = i[0] >> 25;

    auto johni = (i[0] << 8) >> 30;
    auto res = johni * (threehalfs[0] - (x2[0] * johni * johni));
    return res;
}


template <typename T>
BFPDynamic<T> bfp_invsqrtfloat(const BFPDynamic<T> &A) {
    size_t N = A.size();
    // assert(N==B.size());

    vector<double> invsqrtAfloat(A.size());
    for (size_t i=0; i < N; i++){
        int32_t j;
        float x2, y;
        const float threehalfs = 1.5F;
        x2 = A[i] * 0.5F;

        y = A[i];
        j = *(int32_t *) &y;
        j = 0x5f3759df - (j >> 1);
        cout << j << endl;
        y = *(float *) &j;
        y = y * (threehalfs - (x2 * y * y));
        invsqrtAfloat.push_back(y);
    }
    return BFPDynamic<T>(invsqrtAfloat);
}


template <typename T> 
BFPDynamic<T> bfp_log(BFPDynamic<T> &A) {
    BFPDynamic<T> logA;
    size_t N = A.size();

    int min_shifts = numeric_limits<int>::max();

    int bitsdouble = (numeric_limits<T>::digits + 1) * 2;
    int bits = bitsdouble >> 1;

    if (!A.normalized){
        A.normalize();
    }

    typename uTx2<T>::type n;
    typename Tx2<T>::type b;

    for (size_t i=0; i < N; i++){
        n = typename Tx2<T>::type(A[i] >> (A.lazy_list[i] - A.lazy_min)) << bits;
        // if(n <= 8)
        //     return (short)(2 * n);

        // cout << n << endl;
        b = bitsdouble - 1;

        while((b > 2) && (!signbit(n))){
            --b;
            n <<= 1;
        }
        n &= 7 << (bitsdouble - 4);
        n >>= (bitsdouble - 4);

        T n2 = n + (b - 1);
        // cout << int(n2) << endl;
        typename uTx2<T>::type absn2 = n2 ^ (typename Tx2<T>::type(-1 * signbit(n2)));
        // min_shifts = min(bits - floor_log2(absn) - 1, min_shifts);

        min_shifts = min(bits - floor_log2(absn2) - 2, min_shifts);
        logA.lazy_list.push_back(min_shifts);
        logA.push_back((n2 << min_shifts));
    }

    // logA.exponent = A.exponent - min_shifts;
    logA.exponent = (A.exponent >> 1) - min_shifts - (bits >> 1);
    logA.lazy_min = min_shifts;
    logA.normalized = true;
    return logA;
}

short bitlog(unsigned long n)
{
    short b;
    if(n <= 8)
        return (short)(2 * n);

    b=31;
    while((b > 2) && ((long)n > 0)){
        // cout << n << endl;
        --b;
        n <<= 1;
    }
    n &= 0x70000000;
        n >>= 28;
    return (short)n + 8 * (b - 1);
}

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


template <typename T> vector<T> Vlog(const vector<T> &A){
    vector<T> logA(A.size());

    for(int i=0;i<logA.size();i++)
        logA[i] = log2(A[i]);

    return logA;
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
template <typename T> void check_add(T& A, T& B){
    assert(A.size() == B.size());
  size_t N = A.size();
  cout << "******************** Checking addition of: ********************\n"
       << A << " +\n" << B << " = \n\n"
       << A.to_float() << " +\n"  << B.to_float() << "\n\n";

  auto AB = A+B;
  auto AB2 = AB.to_float();
  
  auto ABfloat = A.to_float() + B.to_float();
  cout << "Result of BFP-square root:\n"
       << T(AB) << "\n"
       << T(ABfloat) << " wanted.\n\n";

  cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");


  // auto Afloat = A.to_float(), Bfloat = B.to_float();
  // auto AB = A+B;
  
  // auto ABfloat = A.to_float() + B.to_float();

  // vector<bool> sA(N), sB(N);
  // for(int i=0;i<N;i++){ sA[i] = signbit(A[i]); sB[i] = signbit(B[i]); }

  // cout << "\nFloating point:\n"
  //      << Afloat << " +\n" << Bfloat << " =\n" << ABfloat << "\n\n";

  // cout << "Floating point in BFP-representation:\n"
  //      << T(Afloat) << " +\n" << T(Bfloat) << " =\n" << T(ABfloat) << "\n\n";

  // cout << "Result of BFP-addition:\n"
  //      << T(AB) << "\n"
  //      << T(ABfloat) << " wanted.\n\n";


  // cout << "Result of BFP-addition converted to floating point:\n"
  //      << AB.to_float() << "\n"
  //      << ABfloat << " exact,\n"
  //      << T(ABfloat).to_float() << " wanted.\n\n";

  // cout << "Error compared to exact:\n"
  //      <<  (AB.to_float() - ABfloat) << "\n"
  //      << "Error compared to rounded exact:\n"
  //      << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";

  // cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");

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

template <typename T> void check_mul(T& A, T& B){
  assert(A.size() == B.size());
  size_t N = A.size();

  cout << "******************** Checking multiplication of: ********************\n"
       << A << " *\n" << B << " = \n\n"
       << A.to_float() << " *\n"  << B.to_float() << "\n\n";

  auto Afloat = A.to_float(), Bfloat = B.to_float();
  auto AB = A*B;
  auto ABnorm = AB.to_float();

  auto ABfloat = A.to_float() * B.to_float();
  vector<bool> sA(N), sB(N);
  for(int i=0;i<N;i++){ sA[i] = signbit(A[i]); sB[i] = signbit(B[i]); }

  cout << "\nFloating point:\n"
       << Afloat << " *\n" << Bfloat << " =\n" << ABfloat << "\n\n";

  cout << "Floating point in BFP-representation:\n"
       << T(Afloat) << " *\n" << T(Bfloat) << " =\n" << T(ABfloat) << "\n\n";

  cout << "Result of BFP-multiplication converted to floating point:\n"
       << AB.to_float() << "\n"
       << ABfloat << " exact,\n"
       << T(ABfloat).to_float() << " wanted.\n\n";
    
  cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");
  cout << "Error compared to exact:\n"
       <<  (AB.to_float() - ABfloat) << "\n"
       << "Error compared to rounded exact:\n"
       << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";

  // cout << "AB.exponent: " << AB.exponent << endl;
  // AB.exponent = 0;

  cout << "Result of BFP-multiplication:\n"
       << T(ABnorm) << "\n"
       << T(ABfloat) << " wanted.\n\n";

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


template <typename T> void check_sqrt(T& A){
  size_t N = A.size();
  cout << "******************** Checking square root function of: ********************\n"
       << "sqrt " << A << " = \n\n";

  // auto Afloat = A.to_float();
  auto sqrtA = bfp_sqrt(A);
  cout << sqrtA << endl;
  auto Asqrt = sqrtA.to_float();
  cout << sqrtA << endl;

  auto sqrtAfloat = Vsqrt(A.to_float());
  cout << "Result of BFP-square root:\n"
       << T(Asqrt) << "\n"
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


template <typename T> void check_log(T& A){
  size_t N = A.size();
  cout << "******************** Checking log function of: ********************\n"
       << "sqrt " << A << " = \n\n";

  auto logA = bfp_log(A);
  auto Alog = logA.to_float();
  auto logAfloat = Vlog(A.to_float());
  cout << "Result of BFP-square root:\n"
       << T(Alog) << "\n"
       << T(logAfloat) << " wanted.\n\n";

  cout << "Is the result correct? " << (logA.to_float() == T(logAfloat).to_float()? "Yes.\n" : "No.\n");
}
