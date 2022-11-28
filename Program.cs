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
      RegularMesh rm = new RegularMesh(0, 1, 100);
      ApproximatedSolution aS = new ApproximatedSolution(rm.numberOfSteps, rm.StepLength, new Vector(new double[] { 1, 2 }));
      aS.Solve();
      Console.ReadLine();
    }
  }
}

