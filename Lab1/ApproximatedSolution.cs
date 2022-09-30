using System;
using MathPrimitiveLibrary;

namespace AFLETUNOV_LR1.Lab1
{
  public class ApproximatedSolution
  {
    public Matrix InitMatrix { get; private set; }

    public Matrix ExpandedMatrix { get; private set; }

    public double randomValue { get; }
    public Random r { get; private set; }

    public ApproximatedSolution(int size)
    {
      InitMatrix = new Matrix(size);
      ExpandedMatrix = new Matrix(size + 2, size);
      r = new Random();
      randomValue = (r.NextDouble() - 0.5) * 2;
    }
  }
}
