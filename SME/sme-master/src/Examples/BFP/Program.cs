using System;
using System.IO;
using SME;

namespace BFP
{
	class MainClass
	{
		public static void Main(string[] args)
		{
            new Simulation()
                .BuildCSVFile()
                .BuildGraph()
                .BuildVHDL()
                .BuildCPP()
                .Run(
                    new TopLevelSimulator(),
                    new CarryCalculator(),
                    new AdderWithCarry(),
                    new Pipe()
                );
		}
	}
}