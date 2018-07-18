using SME;
using System;
using System.Threading.Tasks;
using System.Linq;
using SME.VHDL;

// using System.Drawing;
// using System.Drawing.Imaging;



namespace BFP
{

	public class CarryCalculator : Process
	{

		[TopLevelInputBus]
		public interface CarryLine : IBus
		{
			// [InitialValue]
			// bool IsValid { get; set; }

			[FixedArrayLength(100)]
			IFixedArray<byte> Elem { get; set; }
		}


		[InputBus]
        private readonly InputSimulator.InputLine1 DataA = Scope.CreateOrLoadBus<InputSimulator.InputLine1>();
        private readonly InputSimulator.InputLine2 DataB = Scope.CreateOrLoadBus<InputSimulator.InputLine2>();

		[OutputBus]
        private readonly CarryLine DataC = Scope.CreateOrLoadBus<CarryLine>();

		public override async Task Run()
		{
			await ClockAsync();

			// using (var img = System.Drawing.Image.FromFile("data.txt"))
			// using (var bmp = new System.Drawing.Bitmap(img))

			await ClockAsync();
			UInt17 my_int = 16;

			int arrayLength = DataA.Elem[0];
			int arrayLength2 = DataB.Elem[0];
			DataC.Elem[0] = DataA.Elem[0];

			await ClockAsync();
			int size = DataA.Elem[1];
			int size2 = DataB.Elem[1];
			DataC.Elem[1] = DataA.Elem[1];

			await ClockAsync();
			int Aexp = DataA.Elem[2];
			int Bexp = DataB.Elem[2];
			int exp_diff = Aexp - Bexp;
	        int exp_diff2 = size ^ ((exp_diff ^ size) & -Convert.ToInt32(exp_diff < size)); // min(x, y)

			DataC.Elem[2] = (byte) exp_diff2;

			await ClockAsync();

	        bool carry = false;
			for(int i = 3; i < 3 + arrayLength; i++){

				int A = DataA.Elem[i];
				int B = DataB.Elem[i];

				// This can maybe be done only once!

		        long ABi = (((long) A) << exp_diff2) + (long) B;

		        ABi = (ABi >> (exp_diff2)) + ((ABi >> (exp_diff2 - 1)) & 1);
		        int abi = Convert.ToInt32(ABi);
		        carry = (Convert.ToBoolean(A >> (size - 1)) ^ Convert.ToBoolean(abi >> (size - 1))) & (Convert.ToBoolean(B >> (size - 1)) ^ Convert.ToBoolean(abi >> (size - 1)));
		        DataC.Elem[i] = Convert.ToByte(carry);
				await ClockAsync();
			}

			await ClockAsync();

		}
	}
}