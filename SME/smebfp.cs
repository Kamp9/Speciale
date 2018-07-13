using System;
using System.Threading.Tasks;
using System.Collections.Generic;
using SME;


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



    public static BFPStatic plus(BFPStatic A, BFPStatic B) {
        if(A.exp < B.exp) return plus(B, A);
        int size = 32;
        BFPStatic AB = new BFPStatic(A.Len());

        int exp_diff = A.exp - B.exp;
        bool carry = false;

        // int exp_diff2 = min(exp_diff, numeric_limits<T>::digits + 1);
        // int exp_diff2 =  Math.Min(exp_diff, size);
        int exp_diff2 = size ^ ((exp_diff ^ size) & -Convert.ToInt32(exp_diff < size)); // min(x, y)

        // A.size()??
        for(int i=0;i<A.Len();i++){

            long ABi = (((long) A[i]) << exp_diff2) + (long) B[i];

            ABi = (ABi >> (exp_diff2)) + ((ABi >> (exp_diff2 - 1)) & 1);
            int abi = Convert.ToInt32(ABi);
            carry |= (Convert.ToBoolean(A[i] >> (size - 1)) ^ Convert.ToBoolean(abi >> (size - 1))) & (Convert.ToBoolean(B[i] >> (size - 1)) ^ Convert.ToBoolean(abi >> (size - 1)));
        }

        if(carry){
            for(int i=0;i<A.Len();i++){

                bool sign = Convert.ToBoolean(((Convert.ToInt64(A[i]) << exp_diff2) + Convert.ToInt64(B[i])) >> (size - 1));

                long v = Convert.ToInt64(A[i]) << exp_diff2;
                long abs_mask = v >> size - 1;
                long ABi = ((v + abs_mask) ^ abs_mask) + Convert.ToInt64(B[i]);

                // ABi = Math.Abs(ABi >> exp_diff2) - Convert.ToInt64(sign);
                v = ABi >> exp_diff2;
                abs_mask = v >> size - 1;
                ABi = ((v + abs_mask) ^ abs_mask) + Convert.ToInt64(sign);

                bool rounding = Convert.ToBoolean(ABi & 1);
                AB[i] = Convert.ToInt32((ABi >> 1) + Convert.ToInt64(rounding));

            }
        }else{
            for(int i=0;i<A.Len();i++){

                long ABi = (((long) A[i]) << exp_diff2) + (long) B[i];
                bool rounding = Convert.ToBoolean((ABi >> (exp_diff2 - 1)) & 1) && (exp_diff2 > 0);
                bool rounding2 = Convert.ToBoolean((ABi & ((1 << exp_diff2) -1)) == (1 << (exp_diff2 - 1))) && Convert.ToBoolean(ABi >> (size - 1));

                AB[i] = Convert.ToInt32((ABi >> exp_diff2) + Convert.ToInt32(rounding) - Convert.ToInt32(rounding2));

            }
        }

        AB.exp = A.exp + Convert.ToInt32(carry);
        return AB;
    }

[ClockedProcess]
public class BFP : SimpleProcess {
    // public interface IOperand : IBus {
    //     [InitialValue(0)]
    //     int OP { get; set; }
    // }

    public interface IInput : IBus {
        [InitialValue(0)]
        BFPStatic IN { get; set; }
    }

    public interface IResult : IBus {
        BFPStatic OUT { get; set; }
    }

    [InputBus]
    IInput Input1 = Scope.CreateBus<IInput>();
    // ISignal Input2 = Scope.CreateBus<ISignal>();

    [OutputBus]
    IResult Output = Scope.CreateBus<IResult>();

    protected override void OnTick() {
        Output.OUT = plus(Input1.IN, Input1.IN);
    }
}


[ClockedProcess]
public class BFPrunner : Process {

    [InputBus, OutputBus]
    private readonly BFP.IInput Input = Scope.CreateOrLoadBus<BFP.IInput>();

    // private readonly Q1.IOperand Operand = Scope.CreateOrLoadBus<Q1.IOperand>();
    // private readonly Q1.ISignal Signal = Scope.CreateOrLoadBus<Q1.ISignal>();
    // Operand.A = rand;

    [OutputBus]
    private readonly BFP.IResult Result = Scope.CreateOrLoadBus<BFP.IResult>();

    public async override Task Run() {
        await ClockAsync();
        // List<int> a = new List<int>{1,2,3};
        // List<int> b = new List<int>{4,5,6};
        // BFPStatic A = new BFPStatic(a, 2);
        // BFPStatic B = new BFPStatic(b, 2);
        // Input1.IN = 4; //A;
        Result.OUT = plus(Input.IN, Input.IN);

    }
}

public class SMEBFP {
    static void Main(string[] args) {
        Console.WriteLine("ofjgeagodoh");
        new Simulation()
            // .BuildCSVFile()
            // .BuildGraph()
            // .BuildVHDL()
            .Run(
                new BFP(),
                new BFPrunner()
             );
    }
}
}