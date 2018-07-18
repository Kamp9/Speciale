using SME;
using System;
using System.Threading.Tasks;
using System.Linq;

// using System.Drawing;
// using System.Drawing.Imaging;

namespace BFP
{

	public class InputSimulator : Process
	{

		[TopLevelOutputBus]
		public interface InputLine1 : IBus
		{
			[InitialValue]
			bool IsValid { get; set; }

			[FixedArrayLength(100)]
			IFixedArray<byte> Elem { get; set; }
		}


		[TopLevelOutputBus]
		public interface InputLine2 : IBus
		{
				[InitialValue]
				bool IsValid { get; set; }

			[FixedArrayLength(100)]
			IFixedArray<byte> Elem { get; set; }
		}


		[TopLevelOutputBus]
		public interface InputLine3 : IBus
		{
				[InitialValue]
				bool IsValid { get; set; }

			[FixedArrayLength(100)]
			IFixedArray<byte> Elem { get; set; }
		}


		[TopLevelOutputBus]
		public interface InputLine4 : IBus
		{
				[InitialValue]
				bool IsValid { get; set; }

			[FixedArrayLength(100)]
			IFixedArray<byte> Elem { get; set; }
		}

                               		// Num_elems, size, exponent, elemets.....
	    private readonly byte[] InputA = {4, 32, 2, 2, 3, 4, 101};
	    private readonly byte[] InputB = {4, 32, 2, 6, 7, 8, 100};


		[OutputBus]
        private readonly InputLine1 DataA = Scope.CreateOrLoadBus<InputLine1>();
		[OutputBus]
        private readonly InputLine2 DataB = Scope.CreateOrLoadBus<InputLine2>();

		[OutputBus]
        private readonly InputLine3 DataA2 = Scope.CreateOrLoadBus<InputLine3>();
		[OutputBus]
        private readonly InputLine4 DataB2 = Scope.CreateOrLoadBus<InputLine4>();

		public override async Task Run()
		{
			// Make check that InputA[0] and inputB[0] are same length.
			await ClockAsync();

			// swap A and B if A_exp is smaller than B_exp.

			// using (var img = System.Drawing.Image.FromFile("data.txt"))
			// using (var bmp = new System.Drawing.Bitmap(img))
			{
				await ClockAsync();
				for(int i = 0; i < InputA.Length; i++){
					DataA.Elem[i] = InputA[i];
					DataB.Elem[i] = InputB[i];

					DataB2.Elem[i] = InputB[i];
					DataA2.Elem[i] = InputA[i];

					await ClockAsync();
				}
				// Data.IsValid = true;

				// var pixel = bmp.GetPixel(j, i);
				// Data.Color[0] = pixel.R;
				// Data.Color[1] = pixel.G;
				// Data.Color[2] = pixel.B;


				// Data.IsValid = false;
			}

			await ClockAsync();

		}
	}
}