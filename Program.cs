using System;
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
      Console.WriteLine("Task 1");
      Console.WriteLine("Task 2");
      RegularMesh rm = new RegularMesh(0, 1, 500);
      ApproximatedSolution aS = new ApproximatedSolution(rm.numberOfSteps, rm.StepLength, 1e-2, new Vector(new double[] { 1, 2 }));
     // aS.Solve();
      //Console.ReadLine(); 
      Console.WriteLine("Task 3");
      RegularMesh rm1 = new RegularMesh(0, 1, 25);
      RegularMesh rm2 = new RegularMesh(0, 1, 25);
      // НУ
      Func<double, double> phi = (x => 100* x);
      // истинное управление
      Func<double, double> exactManagement = (x => 1000 * x * x * (1 - x));
      Func<double, double> appManagement = (x => 1000 * Math.Sin(10 * x));
      SoutionBuilder sB = new SoutionBuilder(rm1, rm2, 2, 4, phi, exactManagement);
      Matrix trueU = sB.SolveForward(new Vector(new double[] { }), true);
      //trueU[trueU.Rows-1].Show();
      SoutionBuilder sBAp = new SoutionBuilder(rm1, rm2, 2, 4, phi, appManagement);
      Vector appY = sBAp.Solve(trueU[trueU.Rows - 1], 1e-3);
      Console.ReadLine();
    }
  }
}

