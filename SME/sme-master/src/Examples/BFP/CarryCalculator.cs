using SME;
using System;
using System.Threading.Tasks;
using System.Linq;
// using SME.VHDL;

// using System.Drawing;
// using System.Drawing.Imaging;



namespace BFP
{
    [ClockedProcess]
	public class CarryCalculator : SimpleProcess
	{
        [InitializedBus]
		public interface CarryLine : IBus
		{
			// [InitialValue]
			// bool IsValid { get; set; }

            bool Elem { get; set; }
      		}


		[InputBus]
        private readonly InputLine1 DataA = Scope.CreateOrLoadBus<InputLine1>();
        [InputBus]
        private readonly InputLine2 DataB = Scope.CreateOrLoadBus<InputLine2>();
 
        [InputBus]
        private readonly InputLine5 CountLine1 = Scope.CreateOrLoadBus<InputLine5>();

		[OutputBus]
        private readonly CarryLine DataC = Scope.CreateOrLoadBus<CarryLine>();

        //  public override async Task Run()
        int exp_diff2 = 0;
        int size = 0;

        protected override void OnTick()
        {
            //UInt17 my_int = 16;
            uint count = CountLine1.Elem;

            if (count == 0)
            {
                uint arrayLength = DataA.Elem;
                uint arrayLength2 = DataB.Elem;
                DataC.Elem = false;
            }

            if (count == 1)
            {
                size = (int) DataA.Elem;
                uint size2 = DataB.Elem;
                DataC.Elem = false;
            }

            if (count == 2)
            {
                uint Aexp = DataA.Elem;
                uint Bexp = DataB.Elem;
                uint exp_diff = Aexp - Bexp;
                exp_diff2 = (int) (size ^ ((exp_diff ^ size) & -((exp_diff < size) ? 1 : 0)));
                DataC.Elem = false;
            }

            bool carry = false;
            if (count > 2)
            {
                uint A = DataA.Elem;
                uint B = DataB.Elem;

                long ABi = (((long) A) << exp_diff2) + (long) B;

                ABi = (ABi >> (exp_diff2));

                long ABi2 = ((long)(ABi >> (exp_diff2 - 1)));
                long ABi3 = ABi2 & 1;
                ABi = ABi + ABi3;
                int abi = (int)ABi;
                carry = ((((A >> (size - 1))) ^ (abi >> (size - 1))) & (((B >> (size - 1))) ^ ((abi >> (size - 1))))) != 0;

                DataC.Elem = carry;
            }
		}
	}
}