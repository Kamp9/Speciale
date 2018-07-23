using SME;
using System;
using System.Threading.Tasks;
using System.Linq;

// using System.Drawing;
// using System.Drawing.Imaging;

namespace BFP
{

    public class TopLevelSimulator : SimulationProcess
	{

		[TopLevelInputBus]
		public interface InputLine1 : IBus
		{
		//	[InitialValue]
		//	bool IsValid { get; set; }

            [FixedArrayLength(sizeof(int))]
			IFixedArray<byte> Elem { get; set; }
		}


        [TopLevelInputBus]
		public interface InputLine2 : IBus
		{
            //		[InitialValue]
            //	bool IsValid { get; set; }

            [FixedArrayLength(sizeof(int))]
			IFixedArray<byte> Elem { get; set; }
		}


        [TopLevelInputBus]
		public interface InputLine3 : IBus
		{
            //	[InitialValue]
            //	bool IsValid { get; set; }

            [FixedArrayLength(sizeof(int))]
			IFixedArray<byte> Elem { get; set; }
		}


        [TopLevelInputBus]
		public interface InputLine4 : IBus
		{
            //	[InitialValue]
            //	bool IsValid { get; set; }

            [FixedArrayLength(sizeof(int))]
			IFixedArray<byte> Elem { get; set; }
		}


        [TopLevelInputBus]
        public interface InputLine5 : IBus
        {
            [FixedArrayLength(sizeof(int))]
            IFixedArray<byte> Elem { get; set; }
        }

        [TopLevelInputBus]
        public interface InputLine6 : IBus
        {
            [FixedArrayLength(sizeof(int))]
            IFixedArray<byte> Elem { get; set; }
        }

                               		// Num_elems, size, exponent, elemets.....
	    private readonly byte[] InputA = {8, 32, 4, 2, 3, 4, 101, 123, 74};
	    private readonly byte[] InputB = {8, 32, 2, 6, 7, 8, 100, 51, 41};
        //int extra_elements = 3;

		[OutputBus]
        private readonly InputLine1 DataA = Scope.CreateOrLoadBus<InputLine1>();
		[OutputBus]
        private readonly InputLine2 DataB = Scope.CreateOrLoadBus<InputLine2>();

		[OutputBus]
        private readonly InputLine3 DataA2 = Scope.CreateOrLoadBus<InputLine3>();
		[OutputBus]
        private readonly InputLine4 DataB2 = Scope.CreateOrLoadBus<InputLine4>();

        [OutputBus]
        private readonly InputLine5 CountLine1 = Scope.CreateOrLoadBus<InputLine5>();
        [OutputBus]
        private readonly InputLine6 CountLine2 = Scope.CreateOrLoadBus<InputLine6>();

        [InputBus]
        private readonly AdderWithCarry.OutputLine Output = Scope.CreateOrLoadBus<AdderWithCarry.OutputLine>();


        public override async Task Run()
        // protected override async Task OnTickAsync()
		{
            int carry_index = -1;
            byte[] OutputC = new byte[InputA.Length];

			await ClockAsync();

            for (int i = 0; i < InputA.Length; i++)
            {
                DataA.Elem[0] = InputA[i];
                DataB.Elem[0] = InputB[i];
                DataB2.Elem[0] = InputB[i];
                DataA2.Elem[0] = InputA[i];
                CountLine1.Elem[0] = (byte) i;
                CountLine2.Elem[0] = (byte) i;

                await ClockAsync();

                //if (i > extra_elements)
                //{
                OutputC[i] = Output.Elem[0];
                //}
                // await ClockAsync();

                // Data.IsValid = true;
            }

            //await ClockAsync();
            //  DataA.Elem[0] = 0;
            //  DataB.Elem[0] = 0;
            //  DataB2.Elem[0] = 0;
            //  DataA2.Elem[0] = 0;
            //  CountLine1.Elem[0] = (byte) InputA.Length;
            //  CountLine2.Elem[0] = (byte) InputA.Length;
            await ClockAsync();
            carry_index = Output.Elem[0];

            //await ClockAsync();

            Console.WriteLine(carry_index);
            for (int i = 3; i < OutputC.Length; i++){
                Console.WriteLine(OutputC[i]);
                Console.WriteLine(InputA[i] + InputB[i]);
                Console.WriteLine();
            }
        }
	}
}