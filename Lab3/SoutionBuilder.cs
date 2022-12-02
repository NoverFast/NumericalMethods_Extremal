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
    enum ForwardSweepCoefsEnum
    {
      y = 0,
      alpha = 1,
      beta = 2
    }
    public Func<double, double> ManagementFunc { get; set; }
    public Func<double, double> Phi { get; set; }

    private double pMin = 0;
    private double pMax = 4000;

    public double a { get; set; } = 2;
    public double b { get; set; } = 4;
    private double[,] gridMatrix { get; set; }
    private double[,] exactSolution { get; set; }

    public int StepsX { get; set; }
    public int StepsY { get; set; }
    public double T { get; set; }
    public double L { get; set; }

    private List<double> TriDiagonalSolve(double[,] matrix, double[] freeCoefs)
    {
      int n = StepsX; // Размерность матрицы по шагу пространства (numberOfStepsSpace + 1)
      List<double> answerVector = new List<double>();
      List<List<double>> forwardSweepCoefs = new List<List<double>>();
      double y = matrix[0, 0];
      double alpha = -matrix[0, 1] / matrix[0, 0];
      double beta = freeCoefs[0] / matrix[0, 0];
      forwardSweepCoefs.Add(new List<double>() { y, alpha, beta });
      for (int i = 1; i < n; i++)
      {
        y = matrix[i, i] + matrix[i, i - 1] * forwardSweepCoefs[i - 1][(int)ForwardSweepCoefsEnum.alpha];
        alpha = -matrix[i, i + 1] / y;
        beta = (freeCoefs[i] - matrix[i, i - 1] * forwardSweepCoefs[i - 1][(int)ForwardSweepCoefsEnum.beta]) / y;
        forwardSweepCoefs.Add(new List<double>() { y, alpha, beta });
      }
      y = matrix[n, n] + matrix[n, n - 1] * forwardSweepCoefs[n - 1][(int)ForwardSweepCoefsEnum.alpha];
      beta = (freeCoefs[n] - matrix[n, n - 1] * forwardSweepCoefs[n - 1][(int)ForwardSweepCoefsEnum.beta]) / y;
      forwardSweepCoefs.Add(new List<double>() { y, double.NaN, beta });
      answerVector.Add(forwardSweepCoefs[n][(int)ForwardSweepCoefsEnum.beta]);
      int count = 0;
      for (int i = n - 1; i >= 0; i--)
      {
        answerVector.Add(forwardSweepCoefs[i][(int)ForwardSweepCoefsEnum.alpha] * answerVector[count] + forwardSweepCoefs[i][(int)ForwardSweepCoefsEnum.beta]);
        count++;
      }
      answerVector.Reverse();
      return answerVector;
    }

    /*public void ImplicitScheme()
    {
      Array.Clear(gridMatrix);
      Array.Clear(exactSolution);

      // Начальные условия
      for (int i = 0; i <= numberOfStepsSpace; i++)
      {
        gridMatrix[0, i] = phi(i * gridStepSpace);
        exactSolution[0, i] = ExactSolution(0, i * gridStepSpace);
      }

      // Граничные условия
      for (int i = 1; i <= numberOfStepsTime; i++)
      {
        gridMatrix[i, 0] = psi0(i * gridStepTime);
        exactSolution[i, 0] = ExactSolution(i * gridStepTime, 0);
      }

      for (int i = 1; i <= numberOfStepsTime; i++)
      {
        gridMatrix[i, gridMatrix.GetLength(1) - 1] = psi1(i * gridStepTime);
        exactSolution[i, exactSolution.GetLength(1) - 1] = ExactSolution(i * gridStepTime, 1);
      }

      List<double> currentU = new List<double>();
      for (int t = 1; t <= numberOfStepsTime; t++)
      {
        currentU.Clear();
        currentU = TriDiagonalSolve(BuildMatrixImplicit(), BuildCoefsImplicit(t)); // расчитываем температуру в узлах сетки
        for (int x = 0; x < numberOfStepsSpace; x++)
        {
          gridMatrix[t, x] = currentU[x]; // добавляем полученные значения в исходную матрицу-сетку
          exactSolution[t, x] = ExactSolution(t * gridStepTime, x * gridStepSpace);
        }
      }
      Console.WriteLine("Uniform Measure: " + Helper.UniformMeasure(gridMatrix, exactSolution));
    } */

    public SoutionBuilder(double L, double T, int stepsX, int stepsY, Func<double, double> initialManagement)
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
