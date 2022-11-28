using System;
using ExtremalOptimization.Lab2;
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
      ApproximatedSolution aS = new ApproximatedSolution(rm.numberOfSteps, rm.StepLength, new Vector(new double[] { 1, 2 }));
      //      ApproximatedSolution aS = new ApproximatedSolution(rm.numberOfSteps, rm.StepLength, new Vector(new double[] { 1, 2 }));
      aS.Solve();
      Console.ReadLine();
      Console.WriteLine("Task 3");
      double StepsX = 100;
      double StepsY = 100;
      double stepLengthX = 1e-2;
      double stepLengthY = 1e-2;

    }
  }
}

