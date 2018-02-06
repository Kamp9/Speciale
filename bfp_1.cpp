typedef unsigned char      T8;
typedef unsigned short     T16;
typedef unsigned int       T32;
typedef unsigned long      T64;
// typedef unsigned long long T128;


template <class T>
class BfpArray {
public:
	BfpArray(int exp, std::vector<T> elems) : exponent(exp), elements(elems) {}
	void add(BfpArray);
    int getExp() const {return exponent;}
    std::vector<T> getElems() const {return elements;}

private:
	int exponent;
	std::vector<T> elements;
};
