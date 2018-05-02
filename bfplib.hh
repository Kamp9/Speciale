#pragma once

#include <numeric>
#include <inttypes.h>
#include <stdio.h>

using namespace std;

void print_bit_rep(double doubleValue){
	uint8_t *bytePointer = (uint8_t *)&doubleValue;

	for(size_t index = 0; index < sizeof(double); index++)
	{
	    uint8_t byte = bytePointer[index];

	    for(int bit = 0; bit < 8; bit++)
	    {
	        printf("%d", byte&1);
	        byte >>= 1;
	    }
        printf(" ");
	}
	printf("\n");
}

// Langsom O(n)-version af floor(log2(x)) for heltal: mest betydende satte bit
template <typename T> int floor_log2(T value)
{
  int i=0;
  while(value > 0){
    value >>= 1;
    i++;
  }

  return i - 1;
}

// Should be unit tested since negative should get one off
template <typename T> int calc_shifts(int bits, T num){
	int r = 0;
	num = abs(num);
	while (num) // unroll for more speed...
	{
	  r++;
	  num >>= 1;
	}
	// cout << bits - r - 1 << endl;
	return bits - r - 1;// - 1;
}


// uint32_t v; // find the log base 2 of 32-bit v
// int r;      // result goes here

// static const int MultiplyDeBruijnBitPosition[32] = 
// {
//   0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
//   8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
// };

// v |= v >> 1; // first round down to one less than a power of 2 
// v |= v >> 2;
// v |= v >> 4;
// v |= v >> 8;
// v |= v >> 16;

// r = MultiplyDeBruijnBitPosition[(uint32_t)(v * 0x07C4ACDDU) >> 27];

// static const int MultiplyDeBruijnBitPosition2[32] = 
// {
//   0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 
//   31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
// };
// r = MultiplyDeBruijnBitPosition2[(uint32_t)(v * 0x077CB531U) >> 27];


// template <typename T> msb(const &T x);

// template <> msb<uint8_t>(const uint32_t &value)
// {
//   //  const uint8_t table[256] = {0,1,1,2,...};
//   return table[value];    
// }

// template <> msb<uint16_t>(const uint16_t &value)
// {
//   // sådan cirka
//   int v = value > 0xff? (value >> 8) : value;
//   int m = value > 0xff? 16 | 0;
//   return  msb<uint8_t>(v) | m;
// }
// template <> msb<uint32_t>(const uint32_t &value)
// {
//   ///  .... bla bla bla
// }

// template <> msb<uint64_t>(const uint64_t &value)
// {
//     const int tab64[64] = {
//     63,  0, 58,  1, 59, 47, 53,  2,
//     60, 39, 48, 27, 54, 33, 42,  3,
//     61, 51, 37, 40, 49, 18, 28, 20,
//     55, 30, 34, 11, 43, 14, 22,  4,
//     62, 57, 46, 52, 38, 26, 32, 41,
//     50, 36, 17, 19, 29, 10, 13, 21,
//     56, 45, 25, 31, 35, 16,  9, 12,
//     44, 24, 15,  8, 23,  7,  6,  5};

//     value |= value >> 1;
//     value |= value >> 2;
//     value |= value >> 4;
//     value |= value >> 8;
//     value |= value >> 16;
//     value |= value >> 32;
//     return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
// }

// int msb_naive(const uint64_t &x)
// {
//   int i=0;
//   while(i>=0){
//   }  
// }

// msb(x)


