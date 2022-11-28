using MathPrimitivesLibrary;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ExtremalOptimization.Lab3
{
  // Оптимизация функционала методом скорейшего спуска. 
  public class SoutionBuilder
  {
    public Matrix BuildMatrixA(double t)
    {
      return new Matrix(new double[,] 
      { 
        { Math.Cos(t), t },
        { 1.0 / (1.0 + t), Math.Sin(t) }
      });
    }
  }
}
