using SME;
using System;
using System.Threading.Tasks;
using System.Linq;

// using System.Drawing;
// using System.Drawing.Imaging;

namespace BFP
{

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

		public override async Task Run()
		{

			await ClockAsync();

			await ClockAsync();

			int arrayLength = DataA2.Elem[0];
			int arrayLength2 = DataB2.Elem[0];
			int arrayLength3 = CarryLine.Elem[0];

			await ClockAsync();
			int size = DataA2.Elem[1];
			int size2 = DataB2.Elem[1];
			int size3 = CarryLine.Elem[1];	

			await ClockAsync();
			int Aexp = DataA2.Elem[2];
			int Bexp = DataB2.Elem[2];
			int exp_diff2 = CarryLine.Elem[2];
			await ClockAsync();
			// await ClockAsync();

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

	            bool sign = Convert.ToBoolean(((Convert.ToInt64(A) << exp_diff2) + Convert.ToInt64(B)) >> (size - 1));

	            long v = Convert.ToInt64(A) << exp_diff2;
	            long abs_mask = v >> size - 1;

	            long ABi = ((v + abs_mask) ^ abs_mask) + Convert.ToInt64(B);

	            // ABi = Math.Abs(ABi >> exp_diff2) - Convert.ToInt64(sign);
	            v = ABi >> exp_diff2;
	            abs_mask = v >> size - 1;
	            ABi = ((v + abs_mask) ^ abs_mask) + Convert.ToInt64(sign);

	            bool rounding = Convert.ToBoolean(ABi & 1);

	            OutputC[i] = (byte) (Convert.ToInt32((ABi >> c) + Convert.ToInt64(rounding)));
				await ClockAsync();

			}
			// await ClockAsync();
		
			for(int i = 0; i < OutputC.Length; i++){
				Console.WriteLine(OutputC[i]);
			}
			Console.WriteLine(carry_index);
		}
	}
}