#include <iostream>
#include <vector>
#include <math.h>
#include <numeric>

typedef char  T8;
typedef short T16;
typedef int   T32;
typedef long  T64;
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

	std::vector<bool> carryAB;
	for (int i = 0; i < bfpArrayA.elements.size(); i++) {
	 	// carry(A+B) = (sign(A) XOR sign(A+B)) AND (sign(B) XOR sign(A+B))
		carryAB.push_back((std::signbit(bfpArrayA.elements[i]) ^ std::signbit(bfpArrayA.elements[i] + bfpArrayB.elements[i])) &&
		                  (std::signbit(bfpArrayB.elements[i]) ^ std::signbit(bfpArrayA.elements[i] + bfpArrayB.elements[i])));
	}
	int exponentDiff = bfpArrayA.exponent - bfpArrayB.exponent;
	for (int i = 0; i < bfpArrayA.elements.size(); i++) {
		retElements.push_back((carryAB[i] * pow(2, bfpArrayA.exponent)) + bfpArrayA.elements[i] + (bfpArrayB.elements[i] >> exponentDiff));
	}
	// check if we had a carry and therefore needs to move the exponent
	// may be bad cause if statement
	int carry = 0;
	if (accumulate(carryAB.begin(), carryAB.end(), carry)){
		retBfpArray.exponent = bfpArrayA.exponent + 1;
	}
	retBfpArray.elements = retElements;
	return retBfpArray;
}


int main() {
    std::vector<T32> vect1 {1, 2, 3, 4, 5};
    BfpArray<T32> bfp1;
    bfp1.exponent = 2;
    bfp1.elements = vect1;


	std::vector<T32> vect2 {5, 6, 7, 8, 9};
    BfpArray<T32> bfp2;
    bfp2.exponent = 1;
    bfp2.elements = vect2;
    
    BfpArray<T32> bfp3 = add(bfp1, bfp2);
    for (int elem : bfp3.elements) {
        std::cout << elem << std::endl;
    }
}