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

using namespace std;

typedef int8_t   T8;
typedef int16_t  T16;
typedef int32_t  T32;
typedef int64_t  T64;

int log2_64 (uint64_t value) {
    const int tab64[64] = {
    63,  0, 58,  1, 59, 47, 53,  2,
    60, 39, 48, 27, 54, 33, 42,  3,
    61, 51, 37, 40, 49, 18, 28, 20,
    55, 30, 34, 11, 43, 14, 22,  4,
    62, 57, 46, 52, 38, 26, 32, 41,
    50, 36, 17, 19, 29, 10, 13, 21,
    56, 45, 25, 31, 35, 16,  9, 12,
    44, 24, 15,  8, 23,  7,  6,  5};

    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    value |= value >> 32;
    return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}


template <typename T> ostream &operator<<(ostream &s, const vector<T> &xs){
    s << "{"; for(int i=0;i<xs.size();i++) s << xs[i] << (i+1<xs.size()?"," : ""); s << "}";
    return s;
}


template <typename T, size_t N>
struct BFPStatic: public std::array<T,N>{
    int exponent;

    BFPStatic(int exponent=0) : exponent(exponent) {}
    BFPStatic(const std::array<T,N> &A, int exponent) : std::array<T,N>(A), exponent(exponent) {}
    BFPStatic(const std::vector<double> &V) {
        assert(V.size() == N);
    
        std::array<T,N> &A(*this);

        exponent=std::numeric_limits<T>::min();
        for(int i=0;i<N;i++){
            int e_V  = rint(log2(fabs(V[i])));
            exponent = std::max(exponent,e_V);
        }

        exponent -= std::numeric_limits<T>::digits - 1;
        double power = pow(2.0,-exponent);
        for(int i=0;i<N;i++){
            assert(V[i]*power >= std::numeric_limits<T>::min());
            assert(V[i]*power <= std::numeric_limits<T>::max());
            A[i] = floor(V[i]*power); // round(V[i]*power);
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

template <typename T,size_t N>
BFPStatic<T,N> operator+(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
  // Make sure that e_A >= e_B
    if(A.exponent < B.exponent) return B+A;

    BFPStatic<T,N> AB;
    bitset<N> carryAB;
    // int msbAB[N];
    bool carry = false;
    int msb = 1 << (numeric_limits<T>::digits - 1);
    // int carry_mask = 1 << (msb + 1);
    int exp_diff = A.exponent - B.exponent;

    // Compute machine-word addition, carry-bit vector, and whether the carry bit was set.
    // TODO: Needs to be done more carefully for signed values.
    for(size_t i=0;i<N;i++){
        int64_t ABi = A[i] + (B[i]>>exp_diff);
        T       abi = ABi;

        // TODO: Sort out msbAB vs carryAB
        // msbAB[i]   = (ABi & carry_mask)>>1;
        carryAB[i] = (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
        carry     |= carryAB[i];
    }
  
    AB.exponent = A.exponent + carry;

    if(carry)
        for(size_t i=0;i<N;i++){
            AB[i] = ((A[i] + (B[i]>>exp_diff)) | carryAB[i] << msb) >> 1; // +s
        }
    else
        for(size_t i=0;i<N;i++)
            AB[i] = A[i] + (B[i]>>exp_diff); // +s

    return AB;
}


// Should this be implemented like operator+ instead of using it, in order to save a loop over N?
template <typename T,size_t N>
BFPStatic<T,N> operator-(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    BFPStatic<T,N> nA;
    for (size_t i = 0; i < N; i++){ nA[i] = -1 * A[i]; }
    nA.exponent = A.exponent;
    return nA + B;
}


template <typename T,size_t N>
BFPStatic<T,N> operator*(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    // Make sure that e_A >= e_B
    if(A.exponent < B.exponent) return B*A;

    BFPStatic<T,N> AB;
    bitset<N> signAB;
    int msb = numeric_limits<T>::digits;
    uint64_t max = 0;
    int exp_diff = A.exponent - B.exponent;
    for (size_t i = 0; i < N; i++) {
        // Find result signbit for AB
        // A[i] >> (sizeof(T) * 8 - 1 is equivalent to signbit(A[i])
        signAB[i] = (A[i] >> (sizeof(T) * 8 - 1)) ^ (B[i] >> (sizeof(T) * 8 - 1));

        uint64_t ABi = abs(A[i]) * (abs(B[i]) >> exp_diff);
        // ABi - ((ABi - max) & ((ABi - max) >> (sizeof(T) * 8 - 1))) is equivalent to std::max(ABi, max)
        max = std::max(ABi, max);
    }
    int shifts = msb * 2 + 1 - log2_64(max);

    for (size_t i = 0; i < N; i++) {
        AB[i] = ((A[i] * (B[i] >> exp_diff)) >> shifts) | (signAB[i] << (msb+1));
    }

    AB.exponent = A.exponent + B.exponent + shifts;

    return AB;
}

template <typename T,size_t N>
BFPStatic<T,N> operator/(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    // Make sure that e_A >= e_B
    // if(A.exponent < B.exponent) return B*A;

    BFPStatic<T,N> AB;
    bitset<N> signAB;
    int msb = numeric_limits<T>::digits;
    uint64_t max = 0;
    int exp_diff = A.exponent - B.exponent;
    for (size_t i = 0; i < N; i++) {
        // Find result signbit for AB
        // A[i] >> (sizeof(T) * 8 - 1 is equivalent to signbit(A[i])
        signAB[i] = (A[i] >> (sizeof(T) * 8 - 1)) ^ (B[i] >> (sizeof(T) * 8 - 1));

        uint64_t ABi = abs(A[i]) / (abs(B[i]) >> exp_diff);
        // ABi - ((ABi - max) & ((ABi - max) >> (sizeof(T) * 8 - 1))) is equivalent to std::max(ABi, max)
        max = std::max(ABi, max);
    }
    int shifts = msb * 2 + 1 - log2_64(max);

    for (size_t i = 0; i < N; i++) {
        AB[i] = ((A[i] / (B[i] >> exp_diff)) >> shifts) | (signAB[i] << (msb+1));
    }

    AB.exponent = A.exponent - B.exponent + shifts;

    return AB;
}



template <typename T> vector<double> operator+(const vector<T> &A, const vector<T> &B){
    assert(A.size() == B.size());
    vector<double> AB(A.size());

    for(int i=0;i<AB.size();i++){
        AB[i] = A[i] + B[i];
    }
    return AB;
}

template <typename T> vector<double> operator-(const vector<T> &A, const vector<T> &B){
    assert(A.size() == B.size());
    vector<double> AB(A.size());

    for(int i=0;i<AB.size();i++)
        AB[i] = A[i] - B[i];

    return AB;
}

// template <typename T> vector<double> operator*(const vector<T> &A, const vector<T> &B){
//     assert(A.size() == B.size());
//     vector<T> AB(A.size());

//     for(int i=0;i<AB.size();i++)
//         AB[i] = A[i] * B[i];

//     return AB;
// }

// template <typename T> vector<double> operator/(const vector<T> &A, const vector<T> &B){
//     assert(A.size() == B.size());
//     vector<T> AB(A.size());

//     for(int i=0;i<AB.size();i++)
//         AB[i] = A[i] / B[i];

//     return AB;
// }



template <typename T> void check_add(const T& A, const T& B)
{
  assert(A.size() == B.size());
  size_t N = A.size();
  cout << "******************** Checking addition of: ********************\n"
       << A << " +\n" << B << " = \n\n"
       << A.to_float() << " +\n"  << B.to_float() << "\n\n";

  auto Afloat = A.to_float(), Bfloat = B.to_float();
  auto AB = A+B;
  auto ABfloat = A.to_float() + B.to_float();

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
    
  cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");
  cout << "Error compared to exact:\n"
       <<  (AB.to_float() - ABfloat) << "\n"
       << "Error compared to rounded exact:\n"
       << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";

  cout << "signs:\n"
       << sA << "\n"
       << sB << "\n\n";
}



template <typename T> void check_multi(const T& A, const T& B)
{
  assert(A.size() == B.size());
  size_t N = A.size();
  cout << "******************** Checking multiplication of: ********************\n"
       << A << " +\n" << B << " = \n\n"
       << A.to_float() << " +\n"  << B.to_float() << "\n\n";

  auto Afloat = A.to_float(), Bfloat = B.to_float();
  auto AB = A*B;
  auto ABfloat = A.to_float() * B.to_float();

  vector<bool> sA(N), sB(N);
  for(int i=0;i<N;i++){ sA[i] = signbit(A[i]); sB[i] = signbit(B[i]); }

  cout << "\nFloating point:\n"
       << Afloat << " +\n" << Bfloat << " =\n" << ABfloat << "\n\n";

  cout << "Floating point in BFP-representation:\n"
       << T(Afloat) << " +\n" << T(Bfloat) << " =\n" << T(ABfloat) << "\n\n";

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

  cout << "signs:\n"
       << sA << "\n"
       << sB << "\n\n";
}



// int main()
// {
//   // Test af forskellige cases:
//   // 1) A og B er disjunkte:
//   BFPStatic<int8_t,4> Afp{{{-20.25, 4.5, 29.75, 6.75}}};
//   BFPStatic<int8_t,4> Bfp{{{-20.25, 4.5, 29.75, 6.75}}};
  
//   BFPStatic<int8_t,4> A{{3, 18, 119, 27}, 0};
//   BFPStatic<int8_t,4> B{{2, -79, 98, -104}, 0};

//   BFPStatic<int8_t,4> C = A * B;
//   cout << C << endl;
//   // check_add(A,A);
//   // check_add(Afp,Afp);
//   // check_add(A,B);
  
//   return 0;
// }
