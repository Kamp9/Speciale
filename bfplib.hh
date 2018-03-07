#pragma once

#include <inttypes.h>
#include <stdio.h>

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
}

// Langsom O(n)-version af floor(log2(x)) for heltal
template <typename T> int msb(T value)
{
  int i=0;
  while(value > 0){
    value >>= 1;
    i++;
  }
  return i;
}

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


