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

typedef std::int8_t   T8;
typedef std::int16_t  T16;
typedef std::int32_t  T32;
typedef std::int64_t  T64;

using namespace std;

template <typename T> ostream &operator<<(ostream &s, const vector<T> &xs)
{
  s << "{"; for(int i=0;i<xs.size();i++) s << xs[i] << (i+1<xs.size()?"," : ""); s << "}";
  return s;
}


template <typename T, size_t N>
struct BFPStatic: public std::array<T,N>
{
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
    exponent -= std::numeric_limits<T>::digits;

    double power = pow(2.0,-exponent);

    for(int i=0;i<N;i++){
      assert(V[i]*power >= std::numeric_limits<T>::min() && V[i]*power <= std::numeric_limits<T>::max());
      A[i] = round(V[i]*power);
    }
  }

  std::vector<double> to_float() const
  {
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

template <typename T> vector<T> operator+(const vector<T> &A, const vector<T> &B)
{
  assert(A.size() == B.size());
  vector<T> AB(A.size());

  for(int i=0;i<AB.size();i++)
    AB[i] = A[i] + B[i];

  return AB;
}

template <typename T> vector<T> operator-(const vector<T> &A, const vector<T> &B)
{
  assert(A.size() == B.size());
  vector<T> AB(A.size());

  for(int i=0;i<AB.size();i++)
    AB[i] = A[i] - B[i];

  return AB;
}


template <typename T> void check_add(const T& A, const T& B)
{
  cout << "******************** Checking addition of: ********************\n"
       << A << " +\n" << B << " = \n\n"
       << A.to_float() << " +\n"  << B.to_float() << "\n\n";

  auto Afloat = A.to_float(), Bfloat = B.to_float();
  auto AB = A+B;
  auto ABfloat    = A.to_float() + B.to_float();

  vector<bool> sA(10), sB(10);
  for(int i=0;i<10;i++){ sA[i] = signbit(A[i]); sB[i] = signbit(B[i]); }

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

int main()
{
  // Test af forskellige cases:
  // 1) A og B er disjunkte:
  BFPStatic<int8_t,10> Afp{{{-20.25, 4.5, 29.75, 6.75, 20.5, 18.5, 20.25, 0.25, 27., 21.}}};
  BFPStatic<int8_t,10> A{{-81, 18, 119, 27, 82, 74, 81, 1, 108, 84},-2};
  BFPStatic<int8_t,10> B{{-39, -79, 98, -104, 4, 6, 57, 23, 75, 88},-2};


  check_add(A,A);
  check_add(Afp,Afp);
  check_add(A,B);
  
  return 0; 
}
