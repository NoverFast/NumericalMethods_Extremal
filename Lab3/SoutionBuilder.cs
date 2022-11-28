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
    public Func<double, double> ManagementFunc { get; set; }
    public Func<double, double> Phi { get; set; }

    public SoutionBuilder(int stepsX, int stepsY, Func<double, double> initialManagement)
    {

    }

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
