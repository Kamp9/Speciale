#include <boost/utility/binary.hpp>
#include <iostream>
#include <vector>
#include <bitset>


typedef unsigned char      T8;
typedef unsigned short     T16;
typedef unsigned int       T32;
typedef unsigned long      T64;
// typedef unsigned long long T128;




template <typename T>
struct BfpArray {
	int exponent;
	std::vector<T> elements;
};


template <typename T>
BfpArray<T> add(BfpArray<T> bfpArrayA, BfpArray<T> bfpArrayB){
	BfpArray<T> retBfpArray;
	std::vector<T> retElements;



	for (int i = 0; i < bfpArrayA.elements.size(); i++) {
		retElements.push_back(bfpArrayA.elements[i] + bfpArrayB.elements[i]);
	}

	retBfpArray.exponent = bfpArrayA.exponent + 1;
	retBfpArray.elements = retElements;
	return retBfpArray;
}


int main() {
    std::vector<T32> vect1 {0x3, 0x3,0x3,0x3,0x3};
    BfpArray<T32> bfp1;
    bfp1.exponent = 0x6;
    bfp1.elements = vect1;


	std::vector<T32> vect2 {0x6, 0x6,0x6,0x6,0x6};
    BfpArray<T32> bfp2;
    bfp2.exponent = 0x9;
    bfp2.elements = vect2;
    
    BfpArray<T32> bfp3 = add(bfp1, bfp2);
    for (int elem : bfp3.elements) {
        std::cout << elem << std::endl;
     	std::cout << bfp3.exponent << std::endl;
    }
}