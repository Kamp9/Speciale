using SME;
using System;
using System.Threading.Tasks;
using System.Linq;
using SME.VHDL;

namespace BFP
{

    [TopLevelInputBus, InitializedBus]
    public interface InputLine1 : IBus
    {
        //  [InitialValue]
        //  bool IsValid { get; set; }
        uint Elem { get; set; }
    }


    [TopLevelInputBus, InitializedBus]
    public interface InputLine2 : IBus
    {
        uint Elem { get; set; }
    }


    // [TopLevelInputBus, ClockedBus]
    [InitializedBus]
    public interface InputLine3 : IBus
    {
        //  [InitialValue]
        //  bool IsValid { get; set; }
        uint Elem { get; set; }
    }


    // [TopLevelInputBus, ClockedBus]
    [InitializedBus]
    public interface InputLine4 : IBus
    {
        //  [InitialValue]
        //  bool IsValid { get; set; }
        uint Elem { get; set; }
    }


    [TopLevelInputBus, InitializedBus]
    public interface InputLine5 : IBus
    {
        uint Elem { get; set; }
    }

    [InitializedBus]
    public interface InputLine6 : IBus
    {
        uint Elem { get; set; }
    }

    public class TopLevelSimulator : SimulationProcess
	{

        // Num_elems, size, exponent, elemets.....
        private readonly uint[] InputA = { 8, 32, 1, 2, 3, 4, 101, 123, 74, 2147483640};
        private readonly uint[] InputB = { 8, 32, 0, 6, 7, 7, 100, 51, 41, 100};
        //int extra_elements = 3;

		[OutputBus]
        private readonly InputLine1 DataA = Scope.CreateOrLoadBus<InputLine1>();
		[OutputBus]
        private readonly InputLine2 DataB = Scope.CreateOrLoadBus<InputLine2>();

		/*
        [OutputBus]
        private readonly InputLine3 DataA2 = Scope.CreateOrLoadBus<InputLine3>();
		[OutputBus]
        private readonly InputLine4 DataB2 = Scope.CreateOrLoadBus<InputLine4>();
        */

        [OutputBus]
        private readonly InputLine5 CountLine1 = Scope.CreateOrLoadBus<InputLine5>();

        //[OutputBus]
        //private readonly InputLine6 CountLine2 = Scope.CreateOrLoadBus<InputLine6>();

        [InputBus]
        private readonly AdderWithCarry.OutputLine Output = Scope.CreateOrLoadBus<AdderWithCarry.OutputLine>();


        public override async Task Run()
        // protected override async Task OnTickAsync()
		{
            int carry_index = -1;
            uint[] OutputC = new uint[InputA.Length];

			await ClockAsync();

            for (uint i = 0; i < InputA.Length; i++)
            {
                DataA.Elem = InputA[i];
                DataB.Elem = InputB[i];
               // DataB2.Elem = InputB[i];
               // DataA2.Elem = InputA[i];
                CountLine1.Elem = i;
               // CountLine2.Elem = i;


                OutputC[i] = Output.Elem;
                await ClockAsync();
            }

            await ClockAsync();
            carry_index = (int) Output.Elem;

            Console.WriteLine(carry_index);
            Console.WriteLine();

            for (int i = 3; i < OutputC.Length; i++){
                
                Console.WriteLine(OutputC[i]);
                Console.WriteLine(InputA[i] * 2.0 + InputB[i]);
                Console.WriteLine();
            }
        }
	}
}