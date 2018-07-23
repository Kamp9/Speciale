using SME;
using System;
using System.Threading.Tasks;
using System.Linq;

// using System.Drawing;
// using System.Drawing.Imaging;

namespace BFP
{
	public class AdderWithCarry : SimpleProcess
	{

        [TopLevelOutputBus]
        public interface OutputLine : IBus
        {
            [FixedArrayLength(sizeof(int))]
            IFixedArray<byte> Elem { get; set; }
        }

	    // private readonly byte[] OutputC = new byte[7];

		[InputBus]
        private readonly CarryCalculator.CarryLine CarryLine = Scope.CreateOrLoadBus<CarryCalculator.CarryLine>();

        [InputBus]
        private readonly TopLevelSimulator.InputLine3 DataA2 = Scope.CreateOrLoadBus<TopLevelSimulator.InputLine3>();

        [InputBus]
        private readonly TopLevelSimulator.InputLine4 DataB2 = Scope.CreateOrLoadBus<TopLevelSimulator.InputLine4>();

        [InputBus]
        private readonly TopLevelSimulator.InputLine6 CountLine2 = Scope.CreateOrLoadBus<TopLevelSimulator.InputLine6>();

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
            int count = CountLine2.Elem[0];

            if (count == last_count)
            {
                Output.Elem[0] = (byte)carry_index;
            }
            else
            {
                last_count = count;

                if (count == 0)
                {
                    arrayLength = DataA2.Elem[0];
                    _ = DataB2.Elem[0];
                    Output.Elem[0] = 0;
                }

                if (count == 1)
                {
                    size = DataA2.Elem[0];
                    _ = DataB2.Elem[0];
                    _ = CarryLine.Elem[0];
                    Output.Elem[0] = 0;
                }

                if (count == 2)
                {
                    _ = DataA2.Elem[0];
                    _ = DataB2.Elem[0];
                    exp_diff2 = CarryLine.Elem[0];
                    Output.Elem[0] = 0;
                }

                if (count > 2) //-1
                {

                    int A = DataA2.Elem[0];
                    int B = DataB2.Elem[0];
                    int c = CarryLine.Elem[0];

                    if (c == 1 && !carry_found)
                    {
                        carry_index = count;
                        carry_found = true;
                    }

                    bool sign = ((((Int64)A) << exp_diff2) + B) >> (size - 1) != 0;

                    long v = ((Int64)A) << exp_diff2;
                    long abs_mask = v >> (size - 1);

                    long ABi = ((v + abs_mask) ^ abs_mask) + ((Int64)B);

                    v = ABi >> exp_diff2;
                    abs_mask = v >> size - 1;
                    ABi = ((v + abs_mask) ^ abs_mask) + (sign ? 1 : 0);
                    bool rounding = (ABi & 1) > 0;
                    Output.Elem[0] = (byte)((Int32)((ABi >> c) + (rounding ? 1 : 0)));
                }
            }
   		}
	}
}   