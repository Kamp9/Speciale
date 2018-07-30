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
    public class Pipe : SimpleProcess
    {

        [InputBus]
        private readonly InputLine1 DataA = Scope.CreateOrLoadBus<InputLine1>();
        [InputBus]
        private readonly InputLine2 DataB = Scope.CreateOrLoadBus<InputLine2>();

        [InputBus]
        private readonly InputLine5 CountLine1 = Scope.CreateOrLoadBus<InputLine5>();


        [OutputBus]
        private readonly InputLine3 DataA2 = Scope.CreateOrLoadBus<InputLine3>();
        [OutputBus]
        private readonly InputLine4 DataB2 = Scope.CreateOrLoadBus<InputLine4>();

        [OutputBus]
        private readonly InputLine6 CountLine2 = Scope.CreateOrLoadBus<InputLine6>();

        protected override void OnTick()
        {
            DataA2.Elem = DataA.Elem;
            DataB2.Elem = DataB.Elem;
            CountLine2.Elem = CountLine1.Elem;
        }
    }

}
