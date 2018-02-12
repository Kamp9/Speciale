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


template <typename T, size_t N>
struct BFPStatic: public std::array<T,N>
{
  int exponent;

  BFPStatic(int exponent=0) : exponent(exponent) {}
  BFPStatic(const std::array<T,N> &A, int exponent) : std::array<T,N>(A), exponent(exponent) {}
  BFPStatic(const std::array<double,N> &V) {
    std::array<T,N> &A(*this);

    exponent=std::numeric_limits<T>::min();
    for(int i=0;i<N;i++){
      int e_V  = rint(log2(fabs(V[i])));
      exponent = std::max(exponent,e_V);
    }
    exponent -= std::numeric_limits<T>::digits;

    double power = pow(2.0,-exponent);
    std::cerr << power << std::endl;

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
};

using namespace std;

template <typename T> ostream &operator<<(ostream &s, const vector<T> &xs)
{
  s << "{"; for(int i=0;i<xs.size();i++) s << xs[i] << (i+1<xs.size()?"," : ""); s << "}";
  return s;
}


template <typename T,size_t N>
BFPStatic<T,N> add(const BFPStatic<T,N> &A, const BFPStatic<T,N> &B){
  using namespace std;
  
  BFPStatic<T,N>  AB;
  bitset<N>  carryAB;
  bool carry = false;
  
  // Make sure that e_A >= e_B
  if(A.exponent < B.exponent) return add(B,A);
  int exp_diff = A.exponent - B.exponent;

  // Compute machine-word addition, carry-bit vector, and whether the carry bit was set.
  // TODO: Needs to be done more carefully for signed values.
  for(size_t i=0;i<N;i++){
    AB[i]      = A[i] + (B[i]>>exp_diff);  
    carryAB[i] = (signbit(A[i]) ^ signbit(AB[i])) & (signbit(B[i]) ^ signbit(AB[i]));
    carry     |= carryAB[i];

    printf("%d + (%d >> %d) = %d + %d = %d [(%d ^ %d) ^ (%d ^ %d) = %d]\n",
	   A[i],B[i],exp_diff,A[i],(B[i]>>exp_diff),AB[i],
	   int(signbit(A[i])),int(signbit(AB[i])),int(signbit(B[i])),int(signbit(AB[i])),int(carryAB[i]));
  }

  AB.exponent = A.exponent + carry;

  // If any carry bit was set, we need to shift out the lsb, and set the msb to the carry
  // TODO: How does this work for signed ints?
  int msb = numeric_limits<T>::digits - 1;	

  if(carry)
    for(size_t i=0;i<N;i++)
      AB[i] = ((A[i]+B[i]) >> 1) | (carryAB[i] << msb);

  return AB;
}
  
  
int main()
{
  BFPStatic<int8_t,10> Afp{{{20.25, 4.5, 29.75, 6.75, 20.5, 18.5, 20.25, 17.25, 27., 21.}}};
  BFPStatic<int8_t,10> A{{81, 18, 119, 27, 82, 74, 81, 69, 108, 84},-2};
  BFPStatic<int8_t,10> B{{-39, -79, 98, -104, 4, 6, 57, 23, 75, 88},-2};

  cout << vector<int>(A.begin(),A.end())     << "\n";
  cout << vector<int>(Afp.begin(),Afp.end()) << "\n\n";

  cout << A.to_float()   << "\n";
  cout << Afp.to_float() << "\n\n";  

  auto AA = add(A,A);
  cout << "\n"<<vector<int>(AA.begin(),AA.end())  << "\n"
       << AA.to_float() << "\n\n";
  
  return 0; 
}
