using MathPrimitivesLibrary;
using MathPrimitivesLibrary.Types.Meshes;
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
    public Func<double, double> PhiFunc { get; set; }

    private double pMin { get { return 0; } }
    private double pMax { get { return 4000; } }

    public double Alpha { get; set; }
    public double Beta { get; set; }

    public int StepsX { get; set; }
    public int StepsY { get; set; }

    public double StepsXLength { get; set; }
    public double StepsYLength { get; set; }

    public RegularMesh rmX { get; set; }
    public RegularMesh rmY { get; set; }
    /// <summary>
    /// Температура на сетке
    /// </summary>
    public Matrix Grid { get; set; }
    public Matrix ExcatSolution { get; set; }
    /// <summary>
    /// Управление
    /// </summary>
    public Vector StartManagement { get; set; }

    public SoutionBuilder(RegularMesh rmX, RegularMesh rmY, 
      double aKoef, double bKoef, Func<double, double> phi, Func<double, double> management)
    {
      this.rmX = rmX;
      this.rmY = rmY;
      StepsX = rmX.numberOfSteps;
      StepsY = rmY.numberOfSteps;

      StepsXLength = rmX.StepLength;
      StepsYLength = rmY.StepLength;

      Alpha = aKoef;
      Beta = bKoef;

      PhiFunc = phi;
      ManagementFunc = management;
      
      StartManagement = new Vector(StepsY);
      for (int i = 0; i < StepsY; i++)
      {
        StartManagement[i] = ManagementFunc(i * rmY.StepLength);
      }
    }

    private double RectangleSquare(Vector vec, double step)
    {
      double sum = 0;
      for (int i =0; i < vec.Size; i++)
      {
        sum += vec[i] * rmX.StepLength;
      }
      return sum;
    }

    private double Norm(Vector vec, Vector trueY)
    {
      Vector tmp = new Vector(vec.Size);
      for (int i = 0; i < tmp.Size; i++)
      {
        tmp[i] = Math.Pow(vec[i] - trueY[i], 2);
      }

      return RectangleSquare(tmp, rmX.StepLength);
    }

    public Vector Solve(Vector y, double precision)
    {
      Vector currManagement = StartManagement;
      double norm = double.PositiveInfinity;
      int iterations = 0;

      while (norm > precision)
      {
        iterations++;
        Matrix u = SolveForward(currManagement, false);
        norm = Norm(u[u.Rows - 1], y);

        Console.WriteLine($"Current iteration: {iterations}, Norm: {norm}");

        Matrix psi = SolveBackwards(u, y);
        Vector tmpManagement = new Vector(psi.Rows);

        for (int i =0; i < tmpManagement.Size; i++)
        {
          if (psi[i, psi[i].Size - 1] >= 0)
          {
            tmpManagement[i] = pMin;
          }
          else
          {
            tmpManagement[i] = pMax;
          }
        }

        Vector integralPointsU = new Vector(psi[0].Size);
        for (int i =0; i < integralPointsU.Size; i++)
        {
          integralPointsU[i] = Alpha * Alpha * Beta * psi[i][psi[i].Size - 1] * (tmpManagement[i] - currManagement[i]);
        }

        Vector intergalPointsL = new Vector(u[u.Rows - 1].Size);
        Matrix tmpU = SolveForward(tmpManagement, false);
        for (int i = 0; i < integralPointsU.Size; i++)
        {
          intergalPointsL[i] = Math.Pow(tmpU[tmpU.Rows - 1, i] - u[u.Rows - 1, i], 2);
        }
        double alpha = Math.Min(-0.5 * RectangleSquare(integralPointsU, rmY.StepLength) /
          RectangleSquare(intergalPointsL, rmX.StepLength), 1);
        for (int i = 0; i < currManagement.Size; i++)
        {
          currManagement[i] = currManagement[i] + alpha * (tmpManagement[i] - currManagement[i]);
        }
      }
      Matrix answerU = SolveForward(currManagement, false);
      return answerU[answerU.Rows - 1];
    }


    public Matrix SolveForward(Vector currentManagement, bool isEmpty)
    {
      if (isEmpty)
      {
        currentManagement = StartManagement;
      }
      Matrix u = new Matrix(StepsX, StepsY);
      Matrix A = new Matrix(StepsX - 2, StepsY - 2);
      Vector B = new Vector(StepsX - 2);

      for (int i = 0; i < u[0].Size; i++)
      {
        u[0, i] = PhiFunc(rmX.StepLength * i);
      }

      double k1 = Alpha * Alpha / (rmX.StepLength * rmX.StepLength);
      double k2 = -((2.0 * Alpha * Alpha) / (rmX.StepLength * rmX.StepLength) + 1.0 / rmY.StepLength);
      double k3 = k1;
      double k4 = -1.0 / rmY.StepLength;

      for (int i =1; i < StepsY; i++)
      {
        double p1 = 0;
        A[0, 0] = k1 + k2;
        A[0, 1] = k3;
        B[0] = k4 * u[i - 1, 1] + k1 * p1 * rmX.StepLength;

        for (int j = 1; j < A.Rows-1; j++)
        {
          A[j, j - 1] = k1;
          A[j, j] = k2;
          A[j, j] = k3;
          B[j] = k4 * u[i - 1, j + 1];
        }

        A[StepsX - 3, StepsX - 4] = k1;
        A[StepsX - 3, StepsX - 3] = k2 + k3 * (1.0 / (1.0 + Beta * rmX.StepLength));
        B[StepsX - 3] = k4 * u[i - 1, u[i - 1].Size - 2] -
          k3 * ((1.0 / (1.0 + Beta * rmX.StepLength)) * Beta * rmX.StepLength * currentManagement[i]);
        Vector tmpU = MathPrimitivesLibrary.Solvers.ExactSolvers.TridiagonalSolver.Solve(A, B);
        for (int j = 0; j < tmpU.Size; j++)
        {
          u[i,j + 1] = tmpU[j];
        }
        u[i, 0] = u[i, 1] - p1 * rmX.StepLength;
        u[i, u[i].Size - 1] = u[i, u[i].Size - 2] * (1.0 / (1.0 + Beta * rmX.StepLength))
          + Beta * rmX.StepLength * currentManagement[i];
      }

      return u;
    }

    public Matrix SolveBackwards(Matrix u, Vector y)
    {
      Matrix psi = new Matrix(StepsY, StepsX);
      for (int i =0; i < psi[psi.Rows -1].Size; i++)
      {
        psi[psi.Rows-1, i] = 2 * (u[u.Rows -1, i] - y[i]);
      }

      Matrix A = new Matrix(StepsX - 2, StepsY - 2);
      Vector B = new Vector(StepsX - 2);

      double k1 = -Alpha * Alpha / (rmX.StepLength * rmX.StepLength);
      double k2 = -((2.0 * Alpha * Alpha) / (rmX.StepLength * rmX.StepLength) + 1.0 / rmY.StepLength);
      double k3 = k1;
      double k4 = 1.0 / rmY.StepLength;

      for (int i = StepsY - 2; i >= 0; i--)
      {
        double p1 = 0;
        A[0, 0] = k1 + k2;
        A[0, 1] = k3;
        B[0] = k4 * psi[i + 1, 1] + k1 * p1 * rmX.StepLength;

        for (int j = 1; j < A.Rows - 1; j++)
        {
          A[j, j - 1] = k1;
          A[j, j] = k2;
          A[j, j] = k3;
          B[j] = k4 * psi[i + 1, j + 1];
        }

        A[StepsX - 3, StepsX - 4] = k1;
        A[StepsX - 3, StepsX - 3] = k2 + k3 * (1.0 / (1.0 + Beta * rmX.StepLength));
        B[StepsX - 3] = k4 * psi[i + 1, psi[i + 1].Size - 2];

        Vector tmpPsi = MathPrimitivesLibrary.Solvers.ExactSolvers.TridiagonalSolver.Solve(A, B);
        for (int j = 0; j < tmpPsi.Size; j++)
        {
          psi[i, j + 1] = tmpPsi[j];
        }
        psi[i, 0] = psi[i, 1] - p1 * rmX.StepLength;
        psi[i, psi[i].Size - 1] = psi[i, psi[i].Size - 2] * (1.0 / (1.0 + Beta * rmX.StepLength));
      }
      return psi;
    }

    enum ForwardSweepCoefsEnum
    {
      y = 0,
      alpha = 1,
      beta = 2
    }
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

  }
}

 