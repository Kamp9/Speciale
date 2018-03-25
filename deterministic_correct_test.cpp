#include "bfpstatic2.cpp"

template <typename T> void check(const T& A, const vector<double> &B){
    size_t N = A.size();
    assert(A.size() == B.size());
    
    auto Af = A.to_float();
    for(size_t i = 0; i < N; i++){
        if (!(Af[i] == B[i])){            
            cout << "Test Failed" << endl
            << Af[i] << endl
            << "i: " << i << " !="<< endl
            << B[i] << endl;
            return;
        }
    }
    cout << "Test Passed" << endl;
}

int main(){
    BFPStatic<int8_t, 1> A{{1}, 0};

    check(A * A, (A.to_float() * A.to_float()));

    return 0;
}
