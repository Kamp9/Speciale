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

// library functions for faster computation times
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
            int e_V  = ceil(log2(fabs(V[i])));
            exponent = std::max(exponent,e_V);
        }
        exponent -= std::numeric_limits<T>::digits;

        double power = pow(2.0,-exponent);
        for(int i=0;i<N;i++){
            assert(floor(V[i]*power) >= std::numeric_limits<T>::min());
            assert(floor(V[i]*power) <= std::numeric_limits<T>::max());
            A[i] = floor(V[i]*power);
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

// BFPStatic operators +, -, *, and /
template <typename T,size_t N>
BFPStatic<T,N> operator+(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    // Make sure that e_A >= e_B
    if(A.exponent < B.exponent) return B+A;

    BFPStatic<T,N> AB;
    bitset<N> carryAB;
    bool carry = false;
    int msb = 1 << numeric_limits<T>::digits; // should this be larger, e.i more bits?
    int exp_diff = A.exponent - B.exponent;

    // Compute machine-word addition, carry-bit vector, and whether the carry bit was set.
    // TODO: Needs to be done more carefully for signed values.
    for(size_t i=0;i<N;i++){
        int64_t ABi = A[i] + (B[i]>>exp_diff);
        T       abi = ABi;
        carryAB[i]  = (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
        carry      |= carryAB[i];
    }

    AB.exponent = A.exponent + carry;
    if(carry)
        for(size_t i=0;i<N;i++)
            AB[i] = ((A[i] + (B[i]>>exp_diff)) | carryAB[i] << msb) >> 1; // +s
    else
        for(size_t i=0;i<N;i++)
            AB[i] = A[i] + (B[i]>>exp_diff); // +s
    return AB;
}



// Should this be implemented like operator+ instead of using it, in order to save a loop over N?
// Should this be inplace instead? Cannot since BFPStatic is read only
template <typename T,size_t N>
BFPStatic<T,N> operator-(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    BFPStatic<T,N> nB;
    for (size_t i = 0; i < N; i++){ nB[i] = -1 * B[i]; }
    nB.exponent = B.exponent;

    if(A.exponent < B.exponent)
        return A + nB;
    else
        return nB + A;
}


// exponent=std::numeric_limits<T>::min();
// for(int i=0;i<N;i++){
//     int e_V  = ceil(log2(fabs(V[i])));
//     exponent = std::max(exponent,e_V);
// }
// exponent -= std::numeric_limits<T>::digits;

// double power = pow(2.0,-exponent);

template <typename T,size_t N>
BFPStatic<T,N> operator*(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    // Make sure that e_A >= e_B
    if(A.exponent < B.exponent) return B*A;

    BFPStatic<T,N> AB;
    uint64_t exp = 0;
    int exp_diff = A.exponent - B.exponent;
    
    for (size_t i = 0; i < N; i++) {
        uint64_t ABi = (abs(A[i]) << exp_diff) * (abs(B[i]));
        exp = max(exp, ABi);
    }
    // could check that shifts is not 0    
    int shifts = ceil(log2(exp)) - numeric_limits<T>::digits;
    // cout << ceil(log2(exp)) << endl;
    // cout << numeric_limits<T>::digits << endl;
    // cout << exp << endl;
    // cout << shifts << endl;

    for (size_t i = 0; i < N; i++){
        AB[i] = ((A[i] << exp_diff) * B[i]) >> shifts;
    }
    AB.exponent = A.exponent + B.exponent + shifts - exp_diff;

    return AB;
}

template <typename T,size_t N>
BFPStatic<T,N> operator/(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    // Make sure that e_A >= e_B
    if(A.exponent < B.exponent) return B/A;

    BFPStatic<T,N> AB;
    bitset<N> signAB;
    int msb = numeric_limits<T>::digits;
    uint64_t max = 0;
    int exp_diff = A.exponent - B.exponent;
    for (size_t i = 0; i < N; i++) {
        // Find result signbit for AB
        signAB[i] = (A[i] >> signbit(A[i])) ^ (B[i] >> signbit(B[i]));

        uint64_t ABi = abs(A[i]) / (abs(B[i]) >> exp_diff);
        // ABi - ((ABi - max) & ((ABi - max) >> (sizeof(T) * 8 - 1))) is equivalent to std::max(ABi, max)
        max = std::max(ABi, max);
    }

    int shifts = log2_64(max) - std::numeric_limits<T>::digits + 1;
    for (size_t i = 0; i < N; i++) {
        AB[i] = ((A[i] / (B[i] >> exp_diff)) >> shifts) | (signAB[i] << (msb+1));
    }

    AB.exponent = A.exponent - B.exponent + shifts;

    return AB;
}


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
    
  cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");
  cout << "Error compared to exact:\n"
       <<  (AB.to_float() - ABfloat) << "\n"
       << "Error compared to rounded exact:\n"
       << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";

  cout << "signs:\n"
       << sA << "\n"
       << sB << "\n\n";
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
  // cout << "JOHN" << endl;
  // cout << A.to_float() << endl;
  // cout << B.to_float() << endl;
  // cout << AB << endl;
  // cout << T(ABfloat) << endl;
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

  cout << "signs:\n"
       << sA << "\n"
       << sB << "\n\n";
}

template <typename T> void check_div(const T& A, const T& B){
  assert(A.size() == B.size());
  size_t N = A.size();
  cout << "******************** Checking division of: ********************\n"
       << A << " /\n" << B << " = \n\n"
       << A.to_float() << " /\n"  << B.to_float() << "\n\n";

  auto Afloat = A.to_float(), Bfloat = B.to_float();
  auto AB = A/B;
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

  cout << "signs:\n"
       << sA << "\n"
       << sB << "\n\n";
}
