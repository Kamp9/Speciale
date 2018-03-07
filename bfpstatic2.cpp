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

#include "bfplib.hh"

using namespace std;

template <typename T> ostream &operator<<(ostream &s, const vector<T> &xs){
    s << "{"; for(int i=0;i<xs.size();i++) s << xs[i] << (i+1<xs.size()?"," : ""); s << "}";
    return s;
}

// TODO: Derive double-length type automatically or pass as template parameter.
typedef int64_t  Tx2;

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

template <typename T,size_t N>
BFPStatic<T,N> operator+(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
  using namespace std;
  
  BFPStatic<T,N>  AB;
  bitset<N>  carryAB;
  int msbAB[N];
  bool carry = false;
  int msb        = numeric_limits<T>::digits - 1;
  int carry_mask = (1<<(msb+1));

  // Make sure that e_A >= e_B
  if(A.exponent < B.exponent) return B+A;
  int exp_diff = A.exponent - B.exponent;

  // Compute machine-word addition, carry-bit vector, and whether the carry bit was set.
  // TODO: Needs to be done more carefully for signed values.
  for(size_t i=0;i<N;i++){
    int64_t ABi = A[i] + (B[i]>>exp_diff);
    T       abi = ABi;

    // TODO: Sort out msbAB vs carryAB
    msbAB[i]   = (ABi & carry_mask)>>1;
    carryAB[i] = (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
    carry     |= carryAB[i];
  }

  //  cerr << "msbAB:\n" << vector<int>(msbAB,msbAB+N) << "\n";
  
  AB.exponent = A.exponent + carry;

  if(carry)
    for(size_t i=0;i<N;i++){
      int64_t ABi = A[i] + (B[i]>>exp_diff);
      AB[i]       = (ABi>>1) | msbAB[i]; // +s
    }
  else
    for(size_t i=0;i<N;i++)
      AB[i] = A[i] + (B[i]>>exp_diff); // +s

  return AB;
}

// // BFPStatic operators +, -, *, and /
// template <typename T,size_t N>
// BFPStatic<T,N> operator+(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
//     // Make sure that e_A >= e_B
//     if(A.exponent < B.exponent) return B+A;

//     BFPStatic<T,N> AB;
//     bitset<N> carryAB;
//     bitset<N> sAB;
//     bool carry = false;
//     int msb = numeric_limits<T>::digits;
//     int exp_diff = A.exponent - B.exponent;

//     for(size_t i=0;i<N;i++){
//         Tx2 ABi = A[i] + (B[i] >> exp_diff);
//         T   abi = ABi;
//         carryAB[i]  = (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
//         carry      |= carryAB[i];
//         sAB[i]      = signbit(AB[i]);
//     }
//     AB.exponent = A.exponent + carry;
//     if(carry)   
//         for(size_t i=0;i<N;i++){
//             Tx2 ABi = ((carryAB[i] << msb) + ((A[i] + (B[i] >> exp_diff))));
//             AB[i] = ABi >> 1;
//         }
//     else
//         for(size_t i=0;i<N;i++){
//             // cout << "2" << endl;
//             Tx2 ABi = A[i] + (B[i] >> exp_diff) + ((B[i] >> (exp_diff - 1)) & 1);
//             cout << (A[i] + (B[i] >> exp_diff) + ((B[i] >> (exp_diff - 1)) & 1)) << endl;
//             AB[i] = ABi;
//         }
//     return AB;
// }



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

template <typename T,size_t N>
BFPStatic<T,N> operator*(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    BFPStatic<T,N>  AB;
    Tx2 max_value = 0;

    for (size_t i = 0; i < N; i++) {
      Tx2 ABi = Tx2(A[i]) * B[i];
      max_value = max(max_value, std::abs(ABi));
    }

    // could check that shifts is not 0
    int shifts = floor_log2(max_value) - numeric_limits<T>::digits;
    for (size_t i = 0; i < N; i++){
        AB[i] = (Tx2(A[i]) * B[i]) >> shifts;
    }
    AB.exponent = A.exponent + B.exponent + shifts;

    return AB;
}

template <typename T,size_t N>
BFPStatic<T,N> operator/(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    BFPStatic<T,N> AB;
    int shifts = 0;
    __int128_t exp = 0;
    if(A.exponent >= B.exponent){

        int exp_diff = A.exponent - B.exponent;
        for (size_t i = 0; i < N; i++) {
            __int128_t ABi = (__int128_t(abs(A[i])) << (numeric_limits<T>::digits + 1)) / abs(B[i]);
            exp = max(exp, ABi);
        }

        shifts = floor_log2(exp) - numeric_limits<T>::digits;
        for(size_t i = 0; i < N; i++){
            __int128_t ABi = ((__int128_t(A[i]) << (numeric_limits<T>::digits + 1)) / B[i]) >> shifts;
            AB[i] = ABi;
        }
    }else{
        int exp_diff = B.exponent - A.exponent;
        for(size_t i = 0; i < N; i++) {
            __int128_t ABi = (__int128_t(abs(A[i])) << (numeric_limits<T>::digits +1)) / abs(B[i]);
            exp = max(exp, ABi);
        }
        shifts = floor_log2(exp) - numeric_limits<T>::digits;
        for(size_t i = 0; i < N; i++){
            __int128_t ABi = ((__int128_t(A[i]) << (numeric_limits<T>::digits + 1)) / B[i]) >> shifts;
            AB[i] = ABi;
        }
    }
    AB.exponent = A.exponent - B.exponent + shifts - (numeric_limits<T>::digits + 1);    
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
    
  cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");
  // cout << "Error compared to exact:\n"
  //      <<  (AB.to_float() - ABfloat) << "\n"
  //      << "Error compared to rounded exact:\n"
  //      << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";

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
    
  cout << "Is the result correct? " << (AB.to_float() == T(ABfloat).to_float()? "Yes.\n" : "No.\n");
  // cout << "Error compared to exact:\n"
  //      <<  (AB.to_float() - ABfloat) << "\n"
  //      << "Error compared to rounded exact:\n"
  //      << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";

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
  // cout << "Error compared to exact:\n"
  //      <<  (AB.to_float() - ABfloat) << "\n"
  //      << "Error compared to rounded exact:\n"
  //      << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";

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
  // cout << "Error compared to exact:\n"
  //      <<  (AB.to_float() - ABfloat) << "\n"
  //      << "Error compared to rounded exact:\n"
  //      << (AB.to_float() - T(ABfloat).to_float()) << "\n\n";

  // cout << "signs:\n"
  //      << sA << "\n"
  //      << sB << "\n\n";
}
