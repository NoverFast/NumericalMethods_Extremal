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
      ApproximatedSolution aS = new ApproximatedSolution(rm.numberOfSteps, rm.StepLength, 1e-4, new Vector(new double[] { 1, 2 }));
      aS.Solve();
      Console.ReadLine();
      Console.WriteLine("Task 3");
      RegularMesh rm1 = new RegularMesh(0, 2, 100);
      RegularMesh rm2 = new RegularMesh(0, 2, 100);
      // НУ
      Func<double, double> phi = (x => 10 * x * x * x);
      // истинное управление
      Func<double, double> exactManagement = (x => 100 * x * x);
      SoutionBuilder sB = new SoutionBuilder(rm1, rm2, 1, 1, phi, exactManagement);

    }
  }
}

