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
      RegularMesh rMesh = new RegularMesh(0, 1, 1);
      rMesh.PrintMeshProperties();
      Rectangle rectangle = new Rectangle(rMesh, (x) => x*x);
      Console.WriteLine(rectangle.Left());
      Console.WriteLine(rectangle.Right());
      Console.WriteLine(rectangle.Middle());
      Console.ReadLine();
    }
  }
}

