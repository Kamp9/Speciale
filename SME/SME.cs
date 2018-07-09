using System;
using System.Threading.Tasks;
using System.Collections.Generic;
using SME;


// template <typename T, size_t N>
// struct BFPStatic: public std::array{
//     int exponent;
//     BFPStatic(int exponent=0) : exponent(exponent) {}
//     BFPStatic(std::array &A, int exponent) : std::array(A), exponent(exponent) {}
// };

// public class BFPStatic<T, N> {
//     private static T[] A;
//     private static int exponent;

//     public BFPStatic() {
//         T[] A = new T[] {};
//         int exponent = 0;
//     }

//     public BFPStatic(T[] a, int exp) {
//         T[] A = a;
//         int exponent = exp;
//     }

//     public int exp {
//         get { return exponent; }
//         set { exponent = exp; }
//     }

//     public int Length {
//         get { return A.Length; }
//     }

//     public T this[int i] {
//         get { return A[i]; }
//         set { A[i] = value; }
//     }


public class BFPStatic {
    private List<int> elems = new List<int>();
    private int exponent = 0;


    public BFPStatic(int N) {
        elems = new List<int>( new int[N] );

    }

    public BFPStatic(List<int> a, int exp) {
        elems = new List<int>(a);
        exponent = exp;
    }

    public int this[int i] 
    {
       get { return elems[i]; }
       set { elems[i] = value; }
    }

    public int exp {
        get { return exponent; }
        set { exponent = exp; }
    }

    public int Len() { return elems.Count; }

    

    public static BFPStatic operator+ (BFPStatic A, BFPStatic B) {
        if(A.exp < B.exp) return B+A;

        BFPStatic AB = new BFPStatic(A.Len());

        int exp_diff = A.exp - B.exp;
        bool carry = false;

        // int exp_diff2 = min(exp_diff, numeric_limits<T>::digits + 1);
        int exp_diff2 = Math.Min(exp_diff, 32);


        // A.size()??
        for(int i=0;i<A.Len();i++){

            long ABi = (((long) A[i]) << exp_diff2) + (long) B[i];

            ABi = (ABi >> (exp_diff2)) + ((ABi >> (exp_diff2 - 1)) & 1);
            int abi = Convert.ToInt32(ABi);
            carry |= (Convert.ToBoolean(Math.Sign(A[i])) ^ Convert.ToBoolean(Math.Sign(abi))) & (Convert.ToBoolean(Math.Sign(B[i])) ^ Convert.ToBoolean(Math.Sign(abi)));
        }

        // Everytime we shift something negative odd number, we have to add one in order to simulate the positive numbers

        if(carry){
            for(int i=0;i<A.Len();i++){
                bool sign = Convert.ToBoolean(Math.Sign((Convert.ToInt64(A[i]) << exp_diff2) + Convert.ToInt64(B[i])));
                long ABi = Math.Abs((Convert.ToInt64(A[i]) << exp_diff2)) + Convert.ToInt64(B[i]);
                ABi = Math.Abs(ABi >> exp_diff2) - Convert.ToInt64(sign);
                bool rounding = Convert.ToBoolean(ABi & 1);
                AB[i] = Convert.ToInt32((ABi >> 1) + Convert.ToInt64(rounding));

            }
        }else{
            for(int i=0;i<A.Len();i++){

                long ABi = (((long) A[i]) << exp_diff2) + (long) B[i];
                bool rounding = Convert.ToBoolean((ABi >> (exp_diff2 - 1)) & 1) && (exp_diff2 > 0);
                bool rounding2 = Convert.ToBoolean((ABi & ((1 << exp_diff2) -1)) == (1 << (exp_diff2 - 1))) && Convert.ToBoolean(Math.Sign(ABi));

                AB[i] = Convert.ToInt32((ABi >> exp_diff2) + Convert.ToInt32(rounding) - Convert.ToInt32(rounding2));

            }
        }

        AB.exp = A.exp + Convert.ToInt32(carry);
        return AB;
    }
}


// BFPStatic operator+(const BFPStatic<T, N> &A, const BFPStatic<T, N> &B){
//     if(A.exponent < B.exponent) return B+A;
    
//     // size_t N = A.size();

//     BFPStatic<T, N> AB;

//     int exp_diff = A.exponent - B.exponent;
//     bool carry = false;

//     int exp_diff2 = min(exp_diff, numeric_limits<T>::digits + 1);
//     for(size_t i=0;i<N;i++){
//         Tx2 ABi = (Tx2(A[i]) << exp_diff2) + (B[i]);
//         ABi = (ABi >> (exp_diff2)) + ((ABi >> (exp_diff2 - 1)) & 1);
//         T abi  = ABi;
//         carry |= (signbit(A[i]) ^ signbit(abi)) & (signbit(B[i]) ^ signbit(abi));
//     }
//     // Everytime we shift something negative odd number, we have to add one in order to simulate the positive numbers
//     if(carry){
//         for(size_t i=0;i<N;i++){
//             bool sign = signbit((Tx2(A[i]) << exp_diff2) + B[i]);
//             Tx2 ABi = abs((Tx2(A[i]) << exp_diff2) + B[i]);
//             ABi = -sign ^ (ABi >> exp_diff2);
//             bool rounding = (ABi & 1);
//             AB[i] = (ABi >> 1) + rounding;
//         }
//     }else{
//         for(size_t i=0;i<N;i++){
//             Tx2 ABi = (Tx2(A[i]) << exp_diff2) + (B[i]);
//             bool rounding = ((ABi >> (exp_diff2 - 1)) & 1) && (exp_diff2 > 0);
//             bool rounding2 = ((ABi & ((1 << exp_diff2) -1)) == (1 << (exp_diff2 - 1))) && signbit(ABi);
//             AB[i] = (ABi >> exp_diff2) + rounding - rounding2;
//         }
//     }
//     AB.exponent = A.exponent + carry;

//     return AB;
// }


public class Q1 : SimpleProcess {
    public interface IOperand : IBus {
        [InitialValue(0)]
        int A { get; set; }
    }

    public interface IResult : IBus {
        int C { get; set; }
    }

    public interface ISignal : IBus {
        [InitialValue(0)]
        bool FLAG { get; set; }
    }

    [InputBus]
    IOperand Input1 = Scope.CreateBus<IOperand>();
    ISignal Input2 = Scope.CreateBus<ISignal>();

    [OutputBus]
    IResult Output = Scope.CreateBus<IResult>();

    int temp = 0;
    protected override void OnTick() {
        if (Input2.FLAG) {
            Output.C = temp;
            temp = 0;
        } else {
            temp += Input1.A;
        }
    }
}

[ClockedProcess]
public class Q1Tester : Process {
    public void Test(int a, int b) {
        if (a == b){
            Console.WriteLine("Q1 Test: Passed");
        } else {
            Console.WriteLine("Q1 Test: Failed");
        }
    }

    [OutputBus]
    private readonly Q1.IResult Result = Scope.CreateOrLoadBus<Q1.IResult>();

    [InputBus, OutputBus]
    private readonly Q1.IOperand Operand = Scope.CreateOrLoadBus<Q1.IOperand>();
    private readonly Q1.ISignal Signal = Scope.CreateOrLoadBus<Q1.ISignal>();

    public async override Task Run() {
        int sum = 0;
        for(int i = 1; i <= 100; i++){
            Random rnd = new Random();
            int rand = rnd.Next(100);
            if (i % 10 == 0){
                await ClockAsync();
                Signal.FLAG = true;
                Result.C = Operand.A;
                await ClockAsync();

                Test(Result.C, sum);
                sum = 0;
            } else {
                await ClockAsync();
                Operand.A = rand;
                Signal.FLAG = false;
                sum += rand;
            }
        }
        Console.WriteLine();
    }
}


public class Run {
    static void Main(string[] args) {
        List<int> a = new List<int>{1,2,3};
        List<int> b = new List<int>{4,5,6};
        BFPStatic A = new BFPStatic(a, 2);
        BFPStatic B = new BFPStatic(b, 2);
        BFPStatic C = A+B;

        Console.WriteLine(C[0]);
        Console.WriteLine(C[1]);
        Console.WriteLine(C[2]);
        // using(var simulation = new Simulation()) {
        //     simulation.Run(
        //         new Q1(),
        //         new Q1Tester()
        //     );
        // }
    }
}
