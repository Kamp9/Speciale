using SME;
using System;
using System.Threading.Tasks;
using System.Linq;

// using System.Drawing;
// using System.Drawing.Imaging;

namespace BFP
{
    [ClockedProcess]
	public class AdderWithCarry : Process
	{

		// [TopLevelInputBus]
		// public interface CarryLine : IBus
		// {
		// 	// [InitialValue]
		// 	// bool IsValid { get; set; }

		// 	[FixedArrayLength(100)]
		// 	IFixedArray<byte> Elem { get; set; }
		// }

	    private readonly byte[] OutputC = new byte[7];

		[InputBus]
        private readonly CarryCalculator.CarryLine CarryLine = Scope.CreateOrLoadBus<CarryCalculator.CarryLine>();

		private readonly InputSimulator.InputLine3 DataA2 = Scope.CreateOrLoadBus<InputSimulator.InputLine3>();
		private readonly InputSimulator.InputLine4 DataB2 = Scope.CreateOrLoadBus<InputSimulator.InputLine4>();


        // protected override async Task OnTickAsync()
        public override async Task Run()
        {
            await ClockAsync();
            await ClockAsync();
            await ClockAsync();
            await ClockAsync();

			int arrayLength = DataA2.Elem[0];
			int arrayLength2 = DataB2.Elem[0];
			int arrayLength3 = CarryLine.Elem[0];

            Console.WriteLine("ADDER WITH CARRY 1111");

			await ClockAsync();
			int size = DataA2.Elem[1];
			int size2 = DataB2.Elem[1];
			int size3 = CarryLine.Elem[1];	

			await ClockAsync();
			int Aexp = DataA2.Elem[2];
			int Bexp = DataB2.Elem[2];
			int exp_diff2 = CarryLine.Elem[2];
			await ClockAsync();

			bool carry_found = false;
			int carry_index = -1;

			for(int i = 3; i < 3 + arrayLength; i++){
				int A = DataA2.Elem[i];
				int B = DataB2.Elem[i];
				int c =	CarryLine.Elem[i];

				if(c == 1 && !carry_found){
					carry_index = i;
					carry_found = true;
				}

                //bool sign = Convert.ToBoolean(((Convert.ToInt64(A) << exp_diff2) + Convert.ToInt64(B)) >> (size - 1));
                bool sign = false;

                long v = ((Int64) A) << exp_diff2;
	            long abs_mask = v >> size - 1;

                long ABi = ((v + abs_mask) ^ abs_mask) + ((Int64) B);

	            // ABi = Math.Abs(ABi >> exp_diff2) - Convert.ToInt64(sign);
	            v = ABi >> exp_diff2;
	            abs_mask = v >> size - 1;
                ABi = ((v + abs_mask) ^ abs_mask) + (sign ? 1 : 0); // Convert.ToInt64(sign);

                bool rounding = (ABi & 1) > 0;

                OutputC[i] = (byte) ((Int32) ((ABi >> c) + (rounding ? 1 : 0)));
				await ClockAsync();

			}
		
			for(int i = 0; i < OutputC.Length; i++){
				Console.WriteLine(OutputC[i]);
			}
			Console.WriteLine(carry_index);
		}
	}
}   