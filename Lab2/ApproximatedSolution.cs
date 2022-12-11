using System;
using System.IO;
using MathPrimitivesLibrary;

namespace ExtremalOptimization.Lab2
{
  enum Point
  {
    x = 0,
    y = 1
  }
  public class ApproximatedSolution
  {
    public int NumberOfSteps { get; set; }
    public double StepLength { get; set; }
    public Vector FinalPoint { get; set; }
    
    private Matrix u { get; set; }
    private Matrix x { get; set; }
    private Matrix psi { get; set; }
    private Matrix diffJ { get; set; }
    private Matrix uTmp { get; set; }
    private Matrix xTmp { get; set; }
    private double Norm { get; set; }

    // A = (cos(t), t)
    //     (1/(1+t), sin(t))

    // B = (1-e^(-t), 0)
    //     (0, 1 + sin(2t))
    public ApproximatedSolution(int steps, double stepLength, double norm, Vector finalPoint)
    {
      NumberOfSteps = steps;
      StepLength = stepLength;
      FinalPoint = finalPoint;
      Norm = norm;
      // Управление, каждый элемент - некоторая точка плоскости.
      u = FillInitialApproximation();

      x = new Matrix(NumberOfSteps, 2);
      psi = new Matrix(NumberOfSteps, 2);
      diffJ = new Matrix(NumberOfSteps, 2);
      uTmp = new Matrix(NumberOfSteps, 2);
      xTmp = new Matrix(NumberOfSteps, 2);
    }

    private Matrix FillInitialApproximation()
    {
      Matrix u0 = new Matrix(NumberOfSteps, 2);
      for (int i = 0; i < NumberOfSteps; i++)
      {
        double t = i * StepLength;
        u0[i, (int)Point.x] = 10;
        u0[i, (int)Point.y] = 20 ;
      }
      //u0.Show();
      return u0;
    }

    private Matrix ForwardSolve(Matrix m)
    {
      Matrix x = new Matrix(NumberOfSteps, 2);
      x[0, 0] = 0;
      x[1, 0] = 0;
      for (int i =1; i < NumberOfSteps; i++)
      {
        double t = i * StepLength;
        x[i, (int)Point.x] = x[i - 1, (int)Point.x] + StepLength * (Math.Cos(t) * x[i - 1, (int)Point.x] +
          t * x[i - 1, (int)Point.y] + m[i, (int)Point.x] * (1 - Math.Exp(-t)));
        x[i, (int)Point.y] = x[i - 1, (int)Point.y] + StepLength * (Math.Sin(t) * x[i - 1, (int)Point.y] +
           x[i - 1, (int)Point.x] / (t + 1) + m[i, (int)Point.y] * (1 + Math.Sin(2 * t)));
      }
      //x.Show();
      return x;
    }

    private Matrix BackwardSolve()
    {
      Matrix psi = new Matrix(NumberOfSteps, 2);
      psi[NumberOfSteps - 1, (int)Point.x] = 2 * (x[NumberOfSteps - 1, (int)Point.x] - FinalPoint[(int)Point.x]);
      psi[NumberOfSteps - 1, (int)Point.y] = 2 * (x[NumberOfSteps - 1, (int)Point.y] - FinalPoint[(int)Point.y]);
      for (int i = NumberOfSteps - 2; i >= 0; i--)
      {
        double t = i * StepLength;
        psi[i, (int)Point.x] = psi[i + 1, (int)Point.x] - StepLength * (Math.Cos(t) * psi[i + 1, (int)Point.x] +
          + psi[i + 1, (int)Point.y] / (t + 1));
        psi[i, (int)Point.y] = psi[i + 1, (int)Point.y] - StepLength * (Math.Sin(t) * psi[i + 1, (int)Point.y] +
          + t * psi[i + 1, (int)Point.x]);
      }
      return psi;
    }

    private double Trapezoid()
    {
      double sum = 0;
      sum += Math.Pow(diffJ[0, (int)Point.x] / 2, 2) + Math.Pow(diffJ[0, (int)Point.y] / 2, 2) +
        Math.Pow(diffJ[NumberOfSteps - 1, (int)Point.x] / 2, 2) + Math.Pow(diffJ[NumberOfSteps - 1, (int)Point.y] / 2, 2);
      for (int i = 1; i < NumberOfSteps - 1; i++)
      {
        sum += Math.Pow(diffJ[i, (int)Point.x], 2) + Math.Pow(diffJ[i, (int)Point.y], 2);
      }
      return sum * StepLength;
    }

    private void WriteToFile()
    {
      ClearFile("../../Lab2/Results/CoordX.txt");
      ClearFile("../../Lab2/Results/CoordY.txt");
      StreamWriter swX = new StreamWriter("../../Lab2/Results/CoordX.txt");
      StreamWriter swY = new StreamWriter("../../Lab2/Results/CoordY.txt");
      for (int i = 0; i < x.Rows; i++)
      {
        swX.WriteLine(i + " " + x[i, (int)Point.x]);
        swY.WriteLine(i + " " + x[i, (int)Point.y]);
      }
      swX.Close();
      swY.Close();

      ClearFile("../../Lab2/Results/ManX.txt");
      ClearFile("../../Lab2/Results/ManY.txt");
      StreamWriter swManX = new StreamWriter("../../Lab2/Results/ManX.txt");
      StreamWriter swManY = new StreamWriter("../../Lab2/Results/ManY.txt");
      for (int i = 0; i < u.Rows; i++)
      {
        swManX.WriteLine(i + " " + u[i, (int)Point.x]);
        swManY.WriteLine(i + " " + u[i, (int)Point.y]);
      }
      swManX.Close();
      swManY.Close();
    }

    private void ClearFile(string path)
    {
      File.WriteAllText(path, string.Empty);
    }
    public void Solve()
    {
      double norm = Double.MaxValue; 
      while (norm > Norm)
      {
        norm = Math.Sqrt(Math.Pow(x[NumberOfSteps - 1, (int)Point.x] - FinalPoint[(int)Point.x], 2) +
        Math.Pow(x[NumberOfSteps - 1, (int)Point.y] - FinalPoint[(int)Point.y], 2)); 
        Console.WriteLine("Norm: {0}", norm);
        x = ForwardSolve(u);
        psi = BackwardSolve();
        for (int i =0; i < NumberOfSteps; i++)
        {
          double t = i * StepLength;
          diffJ[i, (int)Point.x] = psi[i, (int)Point.x] * (1 - Math.Exp(-t));
          diffJ[i, (int)Point.y] = psi[i, (int)Point.y] * (1 + Math.Sin(2 * t));
          // uJ = u_k - J'[u_k]
          uTmp[i, (int)Point.x] = u[i, (int)Point.x] - diffJ[i, (int)Point.x];
          uTmp[i, (int)Point.y] = u[i, (int)Point.y] - diffJ[i, (int)Point.y];
        }
        //diffJ.Show();
        // Вновь решаем  изначальную задачу, но для другого управления.
        xTmp = ForwardSolve(uTmp);
        //xJ.Show();
        double trapQuadrature = Trapezoid();
        double dNorm = Math.Pow(xTmp[NumberOfSteps - 1, (int)Point.x] - x[NumberOfSteps - 1, (int)Point.x], 2) +
          Math.Pow(xTmp[NumberOfSteps - 1, (int)Point.y] - x[NumberOfSteps - 1, (int)Point.y], 2);
        double alpha = 1.0 / 2.0 * (trapQuadrature / dNorm);
        //Console.WriteLine(alpha);
        for (int i =0; i < NumberOfSteps; i++)
        {
          u[i, (int)Point.x] -= alpha * diffJ[i, (int)Point.x];
          u[i, (int)Point.y] -= alpha * diffJ[i, (int)Point.y];
        }
        //u.Show();
        Console.WriteLine($"Current point: {x[NumberOfSteps - 1, (int)Point.x]}\t\t{x[NumberOfSteps - 1, (int)Point.y]}");
      }
      WriteToFile();
    }
  }
}
