using System;
using System.Diagnostics;
using ExtremalOptimization.Lab2;
using ExtremalOptimization.Lab3;
using MathPrimitivesLibrary;
using MathPrimitivesLibrary.Types.Meshes;

namespace ExtremalOptimization
{
  internal class Program
  {
    static void Main(string[] args)
    {
      Stopwatch sw = new Stopwatch();
      Console.WriteLine("Task 1");
      Console.WriteLine("Task 2");
      RegularMesh rm = new RegularMesh(0, 1, 500);
      sw.Start();
      ApproximatedSolution aS = new ApproximatedSolution(rm.numberOfSteps, rm.StepLength, 1e-2, new Vector(new double[] { 1, 2 }));
      // aS.Solve();
      sw.Stop();
      Console.WriteLine($"Task 2 was completed in {sw.Elapsed} seconds");
      sw.Reset();
      //Console.ReadLine(); 
      Console.WriteLine("Task 3");
      RegularMesh rm1 = new RegularMesh(0, 1, 50);
      RegularMesh rm2 = new RegularMesh(0, 1, 50);
      // НУ
      Func<double, double> phi = (x => 5 * x);
      // истинное управление
      Func<double, double> exactManagement = (x => 25 * x * x );
      Func<double, double> appManagement = (x => 80 * Math.Exp(x));
      sw.Start();
      SoutionBuilder sB = new SoutionBuilder(rm1, rm2, 1, 1, phi, exactManagement);
      Matrix trueU = sB.SolveForward(new Vector(new double[] { }), true);
      SoutionBuilder sBAp = new SoutionBuilder(rm1, rm2, 1, 1, phi, appManagement);
      Matrix appU = sBAp.Solve(trueU[trueU.Rows - 1], 1e-5);
      sw.Stop();
      Console.WriteLine($"Task 3 was completed in {sw.Elapsed} seconds");
      sw.Reset();
      sBAp.WriteManagement(trueU[trueU.Rows - 1], "../../Lab3/Results/trueMan.txt");
      sBAp.WriteManagement(appU[appU.Rows - 1], "../../Lab3/Results/appMan.txt");
      Console.ReadLine();
    }
  }
}

