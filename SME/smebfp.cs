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


    public static void plus(int i, int A, int B, BFP.OutputLine C, int size) {

        C.Elem[i] = (byte) (A + B);

        // if(A[0] < B[0]) return plus(B, A, size);
        // C.Elems[i.count] = A.Elems[i.count] + B.Elems[i.count];

        // C.Elems[0] = (byte) ((byte) A.Elems[0] + (byte) B.Elems[0]);
 
        // int exp_diff2 = min(exp_diff, numeric_limits<T>::digits + 1);
        // int exp_diff2 =  Math.Min(exp_diff, size);
        // int exp_diff2 = size ^ ((exp_diff ^ size) & -Convert.ToInt32(exp_diff < size)); // min(x, y)

        // // A.size()??
        // for(int i=0;i<A.Len();i++){

        //     long ABi = (((long) A[i]) << exp_diff2) + (long) B[i];

        //     ABi = (ABi >> (exp_diff2)) + ((ABi >> (exp_diff2 - 1)) & 1);
        //     int abi = Convert.ToInt32(ABi);
        //     carry |= (Convert.ToBoolean(A[i] >> (size - 1)) ^ Convert.ToBoolean(abi >> (size - 1))) & (Convert.ToBoolean(B[i] >> (size - 1)) ^ Convert.ToBoolean(abi >> (size - 1)));
        // }

        // if(carry){
        //     for(int i=0;i<A.Len();i++){

        //         bool sign = Convert.ToBoolean(((Convert.ToInt64(A[i]) << exp_diff2) + Convert.ToInt64(B[i])) >> (size - 1));

        //         long v = Convert.ToInt64(A[i]) << exp_diff2;
        //         long abs_mask = v >> size - 1;
        //         long ABi = ((v + abs_mask) ^ abs_mask) + Convert.ToInt64(B[i]);

        //         // ABi = Math.Abs(ABi >> exp_diff2) - Convert.ToInt64(sign);
        //         v = ABi >> exp_diff2;
        //         abs_mask = v >> size - 1;
        //         ABi = ((v + abs_mask) ^ abs_mask) + Convert.ToInt64(sign);

        //         bool rounding = Convert.ToBoolean(ABi & 1);
        //         AB[i] = Convert.ToInt32((ABi >> 1) + Convert.ToInt64(rounding));

        //     }
        // }else{
        //     for(int i=0;i<A.Len();i++){

        //         long ABi = (((long) A[i]) << exp_diff2) + (long) B[i];
        //         bool rounding = Convert.ToBoolean((ABi >> (exp_diff2 - 1)) & 1) && (exp_diff2 > 0);
        //         bool rounding2 = Convert.ToBoolean((ABi & ((1 << exp_diff2) -1)) == (1 << (exp_diff2 - 1))) && Convert.ToBoolean(ABi >> (size - 1));

        //         AB[i] = Convert.ToInt32((ABi >> exp_diff2) + Convert.ToInt32(rounding) - Convert.ToInt32(rounding2));

        //     }
        // }

        // AB.exp = A.exp + Convert.ToInt32(carry);
        // return AB;
    }




[ClockedProcess]
public class BFP : SimpleProcess {

    [TopLevelInputBus]
    public interface InputCounter : IBus
    {
        [InitialValue(0)]
        int count { get; set; }

        [InitialValue(true)]
        bool IsValid { get; set; }
    }


    // [TopLevelInputBus]
    // public interface InputLineA : IBus
    // {
    //     [FixedArrayLength(2)]
    //     IFixedArray<byte> Elems { get; set; }
    // }


    // [TopLevelInputBus]
    // public interface InputLineB : IBus
    // {
    //     [FixedArrayLength(2)]
    //     IFixedArray<byte> Elems { get; set; }
    // }


    [OutputBus]
    // public interface OutputLine : IBus
    // {
    //     int Elem { get; set; }
    //     // [FixedArrayLength(2)]
    //     // IFixedArray<byte> Elems { get; set; }
    // }


    [TopLevelInputBus]
    public interface OutputLine : IBus
    {
        [InitialValue]
        bool IsValid { get; set; }

        [FixedArrayLength(4)]
        IFixedArray<byte> Elem { get; set; }
    }

    
    private readonly byte[] InputA = { 100, 2, 3, 0};
    private readonly byte[] InputB = { 4, 5, 6, 0};

    [InputBus, OutputBus]
    InputCounter Counter = Scope.CreateBus<InputCounter>();

    // [InputBus]
    // InputLineA InputA = Scope.CreateBus<InputLineA>();

    // [InputBus]
    // InputLineB InputB = Scope.CreateBus<InputLineB>();

    [OutputBus]
    OutputLine Output = Scope.CreateBus<OutputLine>();

    protected override void OnTick() {

        if (Counter.IsValid){
            Counter.IsValid = false;

            plus(Counter.count, InputA[Counter.count], InputB[Counter.count], Output, 2);
            Counter.count++;

            Counter.IsValid = true;
        }
    }
}

[ClockedProcess]
public class BFPrunner : Process {

    [OutputBus]
    private readonly BFP.OutputLine Result = Scope.CreateOrLoadBus<BFP.OutputLine>();

    // [InputBus, OutputBus]
    // private readonly BFP.InputLineA InputA = Scope.CreateOrLoadBus<BFP.InputLineA>();
    // private readonly BFP.InputLineB InputB = Scope.CreateOrLoadBus<BFP.InputLineB>();


    public async override Task Run() {
        await ClockAsync();
        for(int i = 0; i < 3; i++){
            await ClockAsync();
            Console.WriteLine(Result.Elem[i]);
        }
    }
}

public class SMEBFP {
    static void Main(string[] args) {
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