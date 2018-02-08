#include <iostream>
#include <vector>
#include <math.h>
#include <numeric>
#include <bitset>

// typedef char  T8;
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
		carryAB.push_back((std::signbit(bfpArrayA.elements[i]) ^ std::signbit(bfpArrayA.elements[i] + bfpArrayB.elements[i])) &
		                  (std::signbit(bfpArrayB.elements[i]) ^ std::signbit(bfpArrayA.elements[i] + bfpArrayB.elements[i])));
		std::cout << bfpArrayA.elements[i] + bfpArrayB.elements[i] << std::endl;
	}

	int exponentDiff = bfpArrayA.exponent - bfpArrayB.exponent;
	for (int i = 0; i < bfpArrayA.elements.size(); i++) {
		retElements.push_back((carryAB[i] * pow(2, bfpArrayA.exponent)) + bfpArrayA.elements[i] + (bfpArrayB.elements[i] >> exponentDiff));
	}

	// check if we had a carry and therefore needs to move the exponent
	// may be bad cause if statement
	retBfpArray.exponent = bfpArrayA.exponent;
	int carry = 0;
	if (accumulate(carryAB.begin(), carryAB.end(), carry)){
		std::cout << "VI ER HER!!" << std::endl;
		retBfpArray.exponent++;
		for (int i = 0; i < retElements.size(); i++) {
			retElements[i] = retElements[i] >> 1;
		}
	}

	retBfpArray.elements = retElements;
	return retBfpArray;
}


template <typename T>
void print(BfpArray<T> bfpArrayA) {
	std::cout<< "BfpArray e: " << bfpArrayA.exponent <<"\n";
	for (int i = 0; i < bfpArrayA.elements.size(); i++) {
		std::string binary = std::bitset< sizeof(T)*8 >(bfpArrayA.elements[i]).to_string(); //to binary
		std::string twocomp = std::bitset< sizeof(T)*8 >(~bfpArrayA.elements[i] + 1).to_string();
		// " " << twocomp
    	std::cout<< i << ": " << bfpArrayA.elements[i] << " " << binary << "\n";
	}
}


int main() {
    std::vector<T16> vect1 {32767, 1, 1, 1, 1};
    BfpArray<T16> bfp1;
    bfp1.exponent = 2;
    bfp1.elements = vect1;
	

	std::vector<T16> vect2 {1, 2, 2, 2, 2};
    BfpArray<T16> bfp2;
    bfp2.exponent = 2;
    bfp2.elements = vect2;
    BfpArray<T16> bfp3 = add(bfp1, bfp2);
    // print(bfp1);
    print(bfp1);
    print(bfp2);
    print(bfp3);
    // for (int elem : bfp3.elements) {

    // }
}