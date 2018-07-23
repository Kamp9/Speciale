using SME;
using System;
using System.Threading.Tasks;
using System.Linq;
// using SME.VHDL;

// using System.Drawing;
// using System.Drawing.Imaging;



namespace BFP
{
	public class CarryCalculator : SimpleProcess
	{

		public interface CarryLine : IBus
		{
			// [InitialValue]
			// bool IsValid { get; set; }

            [FixedArrayLength(sizeof(int))]
			IFixedArray<byte> Elem { get; set; }
		}   


		[InputBus]
        private readonly TopLevelSimulator.InputLine1 DataA = Scope.CreateOrLoadBus<TopLevelSimulator.InputLine1>();
        [InputBus]
        private readonly TopLevelSimulator.InputLine2 DataB = Scope.CreateOrLoadBus<TopLevelSimulator.InputLine2>();
 
        [InputBus]
        private readonly TopLevelSimulator.InputLine5 CountLine1 = Scope.CreateOrLoadBus<TopLevelSimulator.InputLine5>();

		[OutputBus]
        private readonly CarryLine DataC = Scope.CreateOrLoadBus<CarryLine>();

        //  public override async Task Run()
        int exp_diff2 = 0;
        int size = 0;

        protected override void OnTick()
        {
            // UInt17 my_int = 16;
            int count = CountLine1.Elem[0];

            if (count == 0)
            {
                int arrayLength = DataA.Elem[0];
                int arrayLength2 = DataB.Elem[0];
                DataC.Elem[0] = 0;
            }

            if (count == 1)
            {
                size = DataA.Elem[0];
                int size2 = DataB.Elem[0];
                DataC.Elem[0] = 0;
            }

            if (count == 2)
            {
                int Aexp = DataA.Elem[0];
                int Bexp = DataB.Elem[0];
                int exp_diff = Aexp - Bexp;
                exp_diff2 = size ^ ((exp_diff ^ size) & -((exp_diff < size) ? 1 : 0));
                DataC.Elem[0] = 0;
            }

            bool carry = false;
            if (count > 2)
            {
                int A = DataA.Elem[0];
                int B = DataB.Elem[0];

                long ABi = (((long) A) << exp_diff2) + (long) B;

                ABi = (ABi >> (exp_diff2));

                long ABi2 = ((long)(ABi >> (exp_diff2 - 1)));
                long ABi3 = ABi2 & 1;
                ABi = ABi + ABi3;
                int abi = (Int32)ABi;
                //carry = (((A >> (size - 1)) != 0) ^ (abi >> (size - 1)) != 0) && (((B >> (size - 1)) != 0) ^ ((abi >> (size - 1)) != 0));
                carry = ((((A >> (size - 1))) ^ (abi >> (size - 1))) & (((B >> (size - 1))) ^ ((abi >> (size - 1))))) != 0;

                DataC.Elem[0] = (byte)(carry ? 1 : 0);
            }
		}
	}
}