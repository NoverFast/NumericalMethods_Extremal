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
    enum ForwardSweepCoefsEnum
    {
      y = 0,
      alpha = 1,
      beta = 2
    }
    public Func<double, double> ManagementFunc { get; set; }
    public Func<double, double> Phi { get; set; }

    private double pMin { get { return 0; } }
    private double pMax { get { return 4000; } }

    public double a { get; set; }
    public double b { get; set; }

    public int StepsX { get; set; }
    public int StepsY { get; set; }

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
    public Vector Management { get; set; }
    public Vector f { get; set; }

    public SoutionBuilder(RegularMesh rmX, RegularMesh rmY, 
      double aKoef, double bKoef, Func<double, double> phi, Func<double, double> startManagement)
    {
      this.rmX = rmX;
      this.rmY = rmY;
      StepsX = rmX.numberOfSteps;
      StepsY = rmY.numberOfSteps;

      a = aKoef;
      b = bKoef;
      //Grid = new Matrix(StepsX + 1, StepsY + 1);
      //ExcatSolution = new Matrix(StepsX + 1, StepsY + 1);
      Phi = phi;
      ManagementFunc = startManagement;
      
      Management = new Vector(StepsX);
      for (int i =0; i < Management.Size; i++)
      {
        Management[i] = ManagementFunc(i * rmX.StepLength);
      }

      //FillInitialConditions();
    }

    /*private void FillInitialConditions()
    {
      for (int i = 0; i < Management.Size; i++)
      {
        Management[i] = ManagementFunc(i * rmY.StepLength);
      }
      for (int i = 0; i < Grid.Coloumns; i++)
      {
        Grid[0, i] = Phi(i * rmX.StepLength);
        //ExcatSolution[0, i] = Management(i * rmX.StepLength);
      }
      for (int i = 0; i < Grid.Rows; i++)
      {
        Grid[i, 0] = 0;
        //ExcatSolution[i, 0] = 0;
      }
    } */

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
      for (int i =0; i < tmp.Size; i++)
      {
        tmp[i] = Math.Pow(vec[i] - trueY[i], 2);
      }

      return RectangleSquare(tmp, rmX.StepLength);
    }

    public Matrix Solve(Vector y, double precision)
    {
      Vector currManagement = Management;
      double norm = double.PositiveInfinity;
      int iterations = 0;

      while (norm > precision)
      {
        iterations++;
        Matrix u = SolveForward(currManagement);
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
          integralPointsU[i] = a * a * b * psi[i][psi[i].Size - 1] * (tmpManagement[i] - currManagement[i]);
        }

        double integral = RectangleSquare(integralPointsU, rmY.StepLength);

        Vector intergalPointsL = new Vector(u[u.Rows - 1].Size);
        Matrix tmpU = SolveForward(tmpManagement);
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
      Matrix answerU = SolveForward(currManagement);
      return answerU;
    }


    public Matrix SolveForward(Vector currentManagement)
    {
      Matrix u = new Matrix(StepsX, StepsY);
      Matrix A = new Matrix(StepsX - 2, StepsY - 2);
      Vector B = new Vector(StepsX - 2);

      for (int i = 0; i < u[0].Size; i++)
      {
        u[0, i] = Phi(rmX.StepLength * i);
      }

      double k1 = a * a / (rmX.StepLength * rmX.StepLength);
      double k2 = -((2.0 * a * a) / (rmX.StepLength * rmX.StepLength) + 1.0 / rmY.StepLength);
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

        A[StepsX - 1, StepsX - 2] = k1;
        A[StepsX - 1, StepsX - 1] = k2 + k3 * (1.0 / (1.0 + b + rmX.StepLength));
        B[StepsX - 1] = k4 * u[i - 1, u[i - 1].Size - 2] -
          k3 * ((1.0 / (1.0 + b * rmX.StepLength)) * b * rmX.StepLength * currentManagement[i]);
        Vector tmpU = MathPrimitivesLibrary.Solvers.ExactSolvers.TridiagonalSolver.Solve(A, B);
        for (int j = 0; j < tmpU.Size; j++)
        {
          u[i,j + 1] = tmpU[j];
        }
        u[i, 0] = u[i, 1] - p1 * rmX.StepLength;
        u[i, u[i].Size - 1] = u[i, u[i].Size - 2] * (1.0 / (1.0 + b * rmX.StepLength))
          + b * rmX.StepLength * currentManagement[i];
      }

      return u;
    }

    public Matrix SolveBackwards(Matrix u, Vector y)
    {
      return new Matrix(1, 1);
    }
  }
}

 