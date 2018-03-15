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

template <int N> ostream &operator<<(ostream &s, const bitset<N> &bits){
  for(int i=0;i<bits.size();i++) s << bits[i];
  return s;
}

// TODO: Derive double-length type automatically or pass as template parameter.
typedef int64_t Tx2;

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
            double intpart;
            double fractpart = modf(V[i]*power, &intpart);
            // Would really like to take care of this in the operators instead
            // But on the other hand, we would need to make calculations to avoid it
            // cout << V[i]*power << endl;
            if(fractpart == -0.5){
                cout << "rounding is happening!" << endl;
                A[i] = round(V[i]*power) + 1;
            }
            else
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

template <typename T,size_t N>
BFPStatic<T,N> operator+(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    // Make sure that e_A >= e_B
    if(A.exponent < B.exponent) return B+A;

    BFPStatic<T,N> AB;
    bitset<N> carryAB;
    bitset<N> sAB;
    bool carry = false;
    int msb = numeric_limits<T>::digits - 1;
    int exp_diff = A.exponent - B.exponent;

    for(size_t i=0;i<N;i++){
        Tx2 ABi = A[i] + (B[i] >> exp_diff);
        T   abi = ABi;
        carry  |= (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
    }

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


// Bug when we try to round -x.5,
// This gives us the -x+1 instead of -x.
template <typename T,size_t N>
BFPStatic<T,N> operator*(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
    BFPStatic<T,N> AB;
    Tx2 max_value = 0;

    for (size_t i = 0; i < N; i++) {
        Tx2 ABi = Tx2(A[i]) * B[i];
        max_value = max(max_value, std::abs(ABi));
    }

    // could check that shifts is not 0
    int shifts = floor_log2(max_value) - numeric_limits<T>::digits;

    for (size_t i = 0; i < N; i++){
        Tx2 ABi = Tx2(A[i]) * B[i];
        bool rounding = (ABi >> (shifts - 1)) & 1;
        AB[i] = (ABi >> shifts) + rounding;
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
    int shifts = floor_log2(max_value) - numeric_limits<T>::digits;
    for (size_t i = 0; i < N; i++){
        Tx2 ABi = (Tx2(A[i]) << (numeric_limits<T>::digits + 1)) / B[i];
        bool rounding = ((ABi >> (shifts - 1)) & 1);
        bool rounding2 = ((ABi & ((1 << shifts) -1)) == (1 << (shifts - 1))) && signbit(ABi);

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
    int shifts = floor_log2(max_value) - numeric_limits<T>::digits;

    for (size_t i = 0; i < N; i++){
        Tx2 Api = exp((p[0] << p.exponent) * log(A[i] << A.exponent));
        Ap[i] = Api >> shifts;
    }

    cout << shifts << endl;
    Ap.exponent = A.exponent + shifts - numeric_limits<T>::digits;
    return Ap;
}

//https://stackoverflow.com/questions/19611198/finding-square-root-without-using-sqrt-function

template <typename T, size_t N >
BFPStatic<T, N> bfp_sqrt(const BFPStatic<T, N> &A) {
    BFPStatic<T,N> sqrtA;
    Tx2 max_value = 0;

    for (size_t i=0; i < N; i++){
        Tx2 sqrtAi = A[i] << (numeric_limits<T>::digits + 1 + (A.exponent & 1)); // << 1
        Tx2 y = (sqrtAi >> 1); // (1 + sqrtAi) / 2;
        Tx2 z = 0;
        Tx2 y1;
        while (y != z) {
            cout << z << endl;
            z = y;
            y1 = (y + sqrtAi / y); // + (sqrtAi & 1);
            y = (y1 >> 1) + (y1 & 1);
        }
        max_value = max(max_value, y - (y1 & 1));
    }

    int shifts = floor_log2(max_value) - numeric_limits<T>::digits;
    cout << shifts << endl;
    for (size_t i=0; i < N; i++){
        Tx2 sqrtAi = A[i] << (numeric_limits<T>::digits + 1 + (A.exponent & 1)); // " "+ 1// << 1
        Tx2 y = (sqrtAi >> 1); // (1 + sqrtAi) / 2;
        Tx2 z = 0;
        // cout << sqrtAi << endl;
        Tx2 y1;
        while (y != z) {
            z = y;
            // cout << y + float(sqrtAi) / y << endl;
            y1 = (y + sqrtAi / y);
            y = (y1 >> 1) + (y1 & 1);
        }
        y -= (y1 & 1);
        // Vi skal have if(shifts) branches udenom for løkken i stedet
        // bool rounding2 = sqrtAi // ((sqrtAi & ((1 << shifts) -1)) == (1 << (shifts - 1)));
        if (shifts > 0){
            sqrtA[i] = (y >> shifts) + ((y >> (shifts - 1) & 1));
        } else {
            sqrtA[i] = y << abs(shifts);
        }
        sqrtA.exponent = ((A.exponent + shifts - numeric_limits<T>::digits) >> 1) - signbit(shifts); // (A.exponent & 1);

    }
    return sqrtA;
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
  auto sqrtA = bfp_sqrt(A);
  auto sqrtAfloat = Vsqrt(A.to_float());

  cout << "Result of BFP-power:\n"
       << sqrtA << "\n"
       << T(sqrtAfloat) << " wanted.\n\n";

  cout << "Is the result correct? " << (sqrtA.to_float() == T(sqrtAfloat).to_float()? "Yes.\n" : "No.\n");
}
