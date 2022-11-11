using ExtremalOptimization.Lab1;
using MathPrimitivesLibrary.Types.Quadratures;
using MathPrimitivesLibrary.Types.Meshes;
using System;

namespace ExtremalOptimization
{
  internal class Program
  {
    static void Main(string[] args)
    {
      //ApproximatedSolution aprx = new ApproximatedSolution(10);
      RegularMesh rMesh = new RegularMesh(0, 5, 100);
      rMesh.PrintMeshProperties();
      Rectangle rectangle = new Rectangle(rMesh, (x) => x*x);
      Simpson simpson = new Simpson(rMesh, (x) => x * x);
      Trapezoid trap = new Trapezoid(rMesh, (x) => x * x);
      Console.WriteLine(rectangle.Left());
      Console.WriteLine(rectangle.Right());
      Console.WriteLine(rectangle.Middle());
      Console.WriteLine(trap.Calculate());
   //   Console.WriteLine(simpson.Calculate());
      //rMesh.PrintMeshData();
      Console.ReadLine();
    }
  }
}

