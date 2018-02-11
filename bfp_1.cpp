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


// Should the expoent be able to decrease?
template <typename T>
BfpArray<T> add(BfpArray<T> bfpArrayA, BfpArray<T> bfpArrayB){
	BfpArray<T> retBfpArray;
    std::vector<T> retElements;
	retBfpArray.exponent = bfpArrayA.exponent;

	// make vector of carrybits from A + B
	std::vector<bool> carryAB;
	for (int i = 0; i < bfpArrayA.elements.size(); i++) {
	 	// carry(A+B) = (sign(A) XOR sign(A+B)) AND (sign(B) XOR sign(A+B))
		carryAB.push_back((std::signbit(bfpArrayA.elements[i]) ^ std::signbit((T) (bfpArrayA.elements[i] + bfpArrayB.elements[i]))) &
		                  (std::signbit(bfpArrayB.elements[i]) ^ std::signbit((T) (bfpArrayA.elements[i] + bfpArrayB.elements[i]))));
	}

	int exponentDiff = bfpArrayA.exponent - bfpArrayB.exponent;
	int carry = 0;

	// check if the sum of carryAB is larger than 0
	if (accumulate(carryAB.begin(), carryAB.end(), carry)) {
		// check if e_A >= e_B
		if (exponentDiff >= 0) {
			retBfpArray.exponent++;
			for (int i = 0; i < bfpArrayA.elements.size(); i++) {
				// std::cout << (carryAB[i] * pow(2, sizeof(T)*8) + bfpArrayA.elements[i] + (bfpArrayB.elements[i] >> exponentDiff)) << std::endl;
				retElements.push_back((bfpArrayA.elements[i] + (bfpArrayB.elements[i] >> exponentDiff)) >> 1);
			}
		} else {
			retBfpArray.exponent = bfpArrayB.exponent + 1;
			for (int i = 0; i < bfpArrayA.elements.size(); i++) {
				// std::cout << (carryAB[i] * pow(2, sizeof(T)*8) + bfpArrayA.elements[i] + (bfpArrayB.elements[i] >> exponentDiff)) << std::endl;
				retElements.push_back((bfpArrayB.elements[i] + (bfpArrayA.elements[i] >> abs(exponentDiff))) >> 1);
			}
		}
	} else {
		// check if e_A >= e_B
		if (exponentDiff >= 0) {
			for (int i = 0; i < bfpArrayA.elements.size(); i++) {
				retElements.push_back(bfpArrayA.elements[i] + (bfpArrayB.elements[i] >> exponentDiff));
			}
		} else {
			for (int i = 0; i < bfpArrayA.elements.size(); i++) {
				retElements.push_back(bfpArrayB.elements[i] + (bfpArrayA.elements[i] >> abs(exponentDiff)));
			}			
		}
	}

	retBfpArray.elements = retElements;
	return retBfpArray;
}


template <typename T>
BfpArray<T> multiply(BfpArray<T> bfpArrayA, BfpArray<T> bfpArrayB){
	BfpArray<T> retBfpArray;
    std::vector<T> retElements;
	retBfpArray.exponent = bfpArrayA.exponent + bfpArrayB.exponent;

	for (int i = 0; i < bfpArrayA.elements.size(); i++) {
		retElements.push_back(bfpArrayA.elements[i] * bfpArrayB.elements[i]);
	}

	for (int i = 0; i < bfpArrayA.elements.size(); i++) {
		retElements.push_back(bfpArrayA.elements[i] * bfpArrayB.elements[i]);
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
    	std::cout<< i << ":\t" << bfpArrayA.elements[i] << "\t" << bfpArrayA.elements[i] * pow(2, bfpArrayA.exponent) << "\t" << binary << "\n";
	}
}


template <typename T>
void test(BfpArray<T> bfpArrayA, BfpArray<T> bfpArrayB, BfpArray<T> bfpArrayC) {
	bool testRes = true;
	std::cout << "Test: " << std::endl;
	for (int i = 0; i < bfpArrayA.elements.size(); i++) {
		if (!(bfpArrayA.elements[i] * pow(2, bfpArrayA.exponent) + 
			  bfpArrayB.elements[i] * pow(2, bfpArrayB.exponent) == bfpArrayC.elements[i] * pow(2, bfpArrayC.exponent))) {
			std::cout << "i: " << i << "\t" << bfpArrayA.elements[i] * pow(2, bfpArrayA.exponent) + bfpArrayB.elements[i] * pow(2, bfpArrayB.exponent) << " not " << bfpArrayC.elements[i] * pow(2, bfpArrayC.exponent) << std::endl;
			testRes = false;
		}
	}
	if (testRes) {
		std::cout << "Test Passed" << std::endl;
	} else {
		std::cout << "Test Failed" << std::endl;
	}
}


int main() {
    std::vector<T16> vect1 {10, 20, -40, 32767, 2};
    BfpArray<T16> bfp1;
    bfp1.exponent = 2;
    bfp1.elements = vect1;

	std::vector<T16> vect2 {2, -30, -50, 2, 4};
    BfpArray<T16> bfp2;
    bfp2.exponent = 2;
    bfp2.elements = vect2;
    BfpArray<T16> bfp3 = add(bfp1, bfp2);
    
    std::cout << bfp3.exponent << std::endl;

    print(bfp1);
    print(bfp2);
    print(bfp3);
    // test(bfp1, bfp2, bfp3);

}