using SME;
using System;
using System.Threading.Tasks;
using System.Linq;
using SME.VHDL;

namespace BFP
{
    [ClockedProcess]
	public class AdderWithCarry : SimpleProcess
	{

        [TopLevelOutputBus, InitializedBus]
        public interface OutputLine : IBus
        {
            uint Elem { get; set; }
        }

	    // private readonly byte[] OutputC = new byte[7];

		[InputBus]
        private readonly CarryCalculator.CarryLine CarryLine = Scope.CreateOrLoadBus<CarryCalculator.CarryLine>();

        [InputBus]
        private readonly InputLine3 DataA2 = Scope.CreateOrLoadBus<InputLine3>();

        [InputBus]
        private readonly InputLine4 DataB2 = Scope.CreateOrLoadBus<InputLine4>();

        [InputBus]
        private readonly InputLine6 CountLine2 = Scope.CreateOrLoadBus<InputLine6>();

        [OutputBus]
        private readonly OutputLine Output = Scope.CreateOrLoadBus<OutputLine>();

        int exp_diff2 = -1;
        int size = -1;
        int arrayLength = -1;
        int carry_index = -1;
        bool carry_found = false;
        int last_count = -1;

        protected override void OnTick()
        {
            int count = (int) CountLine2.Elem;

            if (count == last_count)
            {
                
                Output.Elem = (uint) carry_index;
            }
            else
            {
                last_count = (int) count;

                if (count == 0)
                {
                    arrayLength = (int) DataA2.Elem;
                    _ = DataB2.Elem;
                    Output.Elem = 0;
                }

                if (count == 1)
                {
                    size = (int) DataA2.Elem;
                    _ = DataB2.Elem;
                    _ = CarryLine.Elem;
                    Output.Elem = 0;
                }

                if (count == 2)
                {
                    uint Aexp = DataA2.Elem;
                    uint Bexp = DataB2.Elem;
                    uint exp_diff = Aexp - Bexp;
                    exp_diff2 = (int) (size ^ ((exp_diff ^ size) & -((exp_diff < size) ? 1 : 0)));
                    _ = CarryLine.Elem;
                    Output.Elem = 0;
                }

                if (count > 2) //-1
                {

                    uint A = DataA2.Elem;
                    uint B = DataB2.Elem;
                    bool c = CarryLine.Elem;

                    if (c && !carry_found)
                    {
                        carry_index = count - 3;
                        carry_found = true;
                    }


                    bool sign = (((long)(A) << exp_diff2) +(long)B) >> (2*size - 1) != 0;
 
                    long v = (((long) A) << exp_diff2) + (long)B;
                    long mask = v >> (size - 1);
                    long ABi = (v + mask) ^ mask;
                    ABi = (sign ? -1 : 0) ^ (ABi >> exp_diff2);
 
                    bool rounding = (ABi & 1) != 1;
                    Output.Elem = (uint) ((int) (ABi >> (c ? 1 : 0)) + (rounding ? 1 : 0));

                }
            }
   		}
	}
}   


//long v = ((Int64)A) << exp_diff2;
//long abs_mask = v >> size - 1;

//long ABi = ((v + abs_mask) ^ abs_mask) + ((Int64)B);

// ABi = Math.Abs(ABi >> exp_diff2) - Convert.ToInt64(sign);
// v = ABi >> exp_diff2;
// abs_mask = v >> size - 1;
// ABi = ((v + abs_mask) ^ abs_mask) + (sign? 1 : 0); // Convert.ToInt64(sign);

// bool rounding = (ABi & 1) > 0;

// OutputC[i] = (byte) ((Int32) ((ABi >> c) + (rounding? 1 : 0)));


// bool sign = ((((Int64)A) << exp_diff2) + B) >> (size - 1) != 0;
// long v = (((Int64)A) << exp_diff2) + B;
// long abs_mask = v >> (size - 1);

// long ABi = ((v + abs_mask) ^ abs_mask);

// v = ABi >> exp_diff2;
//                    abs_mask = v >> (size - 1);
//                    ABi = ((v + abs_mask) ^ abs_mask) + (sign? 1 : 0);

//                    bool rounding = (ABi & 1) != 0;
//Output.Elem[0] = ((int)(ABi >> (c? 1 : 0) + (rounding? 1 : 0)));
