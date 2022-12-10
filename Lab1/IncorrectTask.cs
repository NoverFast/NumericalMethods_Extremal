using System;
using System.Runtime.ExceptionServices;
using MathPrimitivesLibrary;
using Microsoft.Win32;

namespace ExtremalOptimization.Lab1
{
  public class IncorrectTask
  {
    // Исходные данные. Генерируются случайным образом (по закону равномерного распределения на отрезке [-1, 1])
    public Matrix SourceMatrix { get; private set; }
    public Vector ExactSolution { get; private set; }
    public Vector FreeCoefs { get; private set; }

    // Данные расширенной матрицы и столбца свободных членов
    public Matrix ExpandedMatrix { get; private set; }
    public Vector ExpandedCoefs { get; private set; }

    // Данные расширенной матрицы и столбца свободных членов с добавленным шумом
    public Matrix ApproximatedMatrix { get; private set; }
    public Vector ApproximatedCoefs { get; private set; }

    public Random r { get; private set; }

    public double H { get; set;}
    public double Sigma { get; set; }
    public double LambdaDelta { get; set; }

    // в качестве начального приближения метода ньютона для решения СНАУ можно использовать точное решение
    public IncorrectTask(int size, int seed)
    {
      InitializeData(size, seed); // Инициализируем все начальные данные, с которыми впоследствие будем работать
    }

    private void InitializeData(int size, int seed)
    {
      r = new Random(seed);

      SourceMatrix = GenerateRandomMatrix(size);
      ExactSolution = GenerateRandomSolution(size);
      FreeCoefs = SourceMatrix * ExactSolution;
      ExpandedMatrix = GenerateExpandedMatrix(SourceMatrix);
      ExpandedCoefs = GenerateExpandedCoefs(FreeCoefs);
      ApproximatedMatrix = AddNoise(ExpandedMatrix);
      ApproximatedCoefs = AddNoise(ExpandedCoefs);

#if DEBUG
      Console.WriteLine("Исходная матрица и коэффициенты");
      SourceMatrix.Show();
      FreeCoefs.Show();

      Console.WriteLine("\nРасширенная матрица и коэффициенты без шума");
      ExpandedMatrix.Show();
      ExpandedCoefs.Show();

      Console.WriteLine("\nРасширенная матрица и коэффициенты с шумом");
      ApproximatedMatrix.Show();
      ApproximatedCoefs.Show();
#endif
      //test
    }

    double LFunc()
    {
      double buff0, buff1 = 0, buff2 = 0;
      for (int i =0; i <  ExpandedCoefs.Size; i++)
      {
        buff0 = 0;
        for (int j = 0; j < FreeCoefs.Size; j++)
        {
          buff0 += ApproximatedMatrix[i, j] * ExactSolution[j]; 
        }
        buff0 -= ApproximatedCoefs[i];
        buff1 += Math.Pow(buff0, 2);
      }
      for (int i =0; i < FreeCoefs.Size; i++)
      {
        buff2 += Math.Pow(ExactSolution[i], 2);
      }
      return Math.Sqrt(buff1) + H * Math.Sqrt(buff2) + Sigma;
    }

    private double dL(int k)
    {
      double buff0, buff1 = 0, buff2 = 0, buff3 = 0;
      for (int i = 0; i < ExpandedCoefs.Size; i++)
      {
        buff0 = 0;
        for (int j =0; j < FreeCoefs.Size; j++)
        {
          buff0 += ApproximatedMatrix[i, j] * ExactSolution[j];
        }
        buff0 -= ApproximatedCoefs[i];
        buff1 += ApproximatedMatrix[i, k] * buff0;
        buff2 += Math.Pow(buff0, 2);
      }
      for (int i =0; i < FreeCoefs.Size; i++)
      {
        buff3 += Math.Pow(ExactSolution[i], 2);
      }
      return buff1 / Math.Sqrt(buff2) + H * ExactSolution[k] / Math.Sqrt(buff3);
    }

    private Vector gradL()
    {
      Vector gradient = new Vector(ExactSolution.Size);
      for (int i =0; i < gradient.Size; i++)
      {
        gradient[i] = dL(i);
      }
      return gradient; 
    }

    private Vector CGM(Matrix A, Vector b, double precision)
    {
      Vector xPrev = new Vector(ExpandedMatrix.Rows);
      Vector x = xPrev;
      Vector rPrev = b - A * xPrev;
      Vector r = rPrev;
      Vector zPrev = rPrev;
      Vector z = zPrev;
      double error;
      do
      {
        error = 0;
        double alphaK = rPrev.DotProduct(rPrev) / (ApproximatedMatrix * zPrev).DotProduct(zPrev);
        for (int i = 0; i < A.Rows; i++)
        {
          x[i] = xPrev[i] + alphaK * zPrev[i];
          r[i] = rPrev[i] - alphaK * (A * zPrev)[i];
        }
        double betaK = r.DotProduct(r) / rPrev.DotProduct(rPrev);
        for (int i = 0; i < A.Rows; i++)
        {
          z[i] = r[i] + betaK * zPrev[i];
        }
        error = (x - xPrev).Norm();
        for (int i = 0; i < A.Rows; i++)
        {
          xPrev[i] = x[i];
          rPrev[i] = r[i];
          zPrev[i] = z[i];
        }
      } while (error > precision);
      return x;
    }

    private double FindRho(Vector v)
    {
      double buff0, buff1 = 0, buff2 = 0;
      for (int i =0; i < ApproximatedMatrix.Rows; i++)
      {
        buff0 = 0;
        for (int j =0; j < SourceMatrix.Rows; i++)
        {
          buff0 += ApproximatedMatrix[i, j] * v[j];
        }
        buff0 -= ApproximatedCoefs[i];
        buff1 += Math.Pow(buff0, 2);
      }
      for (int i =0; i < SourceMatrix.Rows; i++)
      {
        buff2 += Math.Pow(v[i], 2);
      }
      return buff1 - Math.Pow(LambdaDelta + Sigma + H * Math.Sqrt(buff2), 2); 
    }

    private void SigmaAndH()
    {
      Matrix aTmp = SourceMatrix - ApproximatedMatrix;
      Vector fTmp = FreeCoefs - ApproximatedCoefs;
      H = 0;
      Sigma = fTmp.Norm();
      Console.WriteLine($"H: {H}, | Sigma: {Sigma}");

    }

    private Vector GradientDescent(Func<Vector, double> f, Func<Vector, Vector> gradF, Vector x0, double precision, int M )
    {
      int k = 0;
      Vector xPrev = x0;
      Vector x = xPrev;
      double tK = 3;
      while (true)
      {
        k += 1;
        Vector dd = gradF(xPrev) * tK;
        Vector dd2 = xPrev - dd;
        x = dd2;
        if (f(x) - f(xPrev) < 0)
        {
          Vector tmp = x - xPrev;
          if (tmp.Norm() <= precision || k + 1 > M)
          {
            return x;
          }
        }
        else
        {
          tK /= 2;
        }
        for (int i =0; i < SourceMatrix.Rows; i++)
        {
          xPrev[i] = x[i];
        }
      }
    }

    private Matrix AddNoise(Matrix m)
    {
      Matrix noiseMatrix = new Matrix(m);
      for (int i =0; i < noiseMatrix.Rows; i++)
      {
        for (int j =0; j < noiseMatrix.Coloumns; j++)
        {
          noiseMatrix[i, j] *= (1 + ((r.NextDouble() - 0.5) * 2) * 0.001);
        }
      }
      return noiseMatrix;
    }

    private void FindLambdaDelta()
    {
      Vector x0 = new Vector(ExactSolution);
      Vector LMin = GradientDescent(LFunc, gradL, x0, 1e-6, 1000);
      LambdaDelta = LFunc(LMin);
      LMin.Show();
      Console.Write($"\nLambda Delta: {LambdaDelta}");
    }
    private void FindSolution()
    {
      double alpha = 1;
      Matrix matrix = ApproximatedMatrix.Transpose() * ApproximatedMatrix + SourceMatrix.IdentityMatrix * alpha;
      Vector vector = ApproximatedMatrix.Transpose() * ApproximatedCoefs;
      Vector answer = CGM(matrix, vector, 1e-6);
      matrix.Show();
      vector.Show();
      // вывод в файл
      //
      double val = FindRho(answer);

      double needAlpha = 0;
      double needValue = 0;
      Vector needAnswer = new Vector(answer.Size);
      while (alpha >= -0.01)
      {
        alpha -= 0.0001;
        matrix = ApproximatedMatrix.Transpose() * ApproximatedMatrix + SourceMatrix.IdentityMatrix * alpha;
        vector = ApproximatedMatrix.Transpose() * ApproximatedCoefs;
        answer = CGM(matrix, vector, 1e-6);
        val = FindRho(answer);
        if (alpha <= 0.25)
        {
          //запись в файл
        }
        if (Math.Abs(val) < 1e-4 && alpha >= 0)
        {
          needAlpha = alpha;
          needValue = val;
          needAnswer = answer;
        }
      }
      // закрыть файл
      Console.WriteLine($"Value: {needValue}");
      Console.WriteLine($"Alpha: {needAlpha}");
      needAnswer.Show();
      Console.WriteLine($"\nAnswer Norm: {needAnswer.Norm()}");
    }
    private Vector AddNoise(Vector v)
    {
      Vector noiseVector = new Vector(v);
      for (int i = 0; i < noiseVector.Size; i++)
      {
        noiseVector[i] *= (1 + ((r.NextDouble() - 0.5) * 2) * 0.01);
      }
      return noiseVector;
    }
    private Matrix GenerateRandomMatrix(int size)
    {
      Matrix m = new Matrix(size);
      double randomValue = 0;
      for (int i =0; i < m.Rows; i++)
      {
        for (int j =0; j < m.Coloumns; j++)
        {
          randomValue = (r.NextDouble() - 0.5) * 2;
          m[i, j] = randomValue;
        }
      }
      return m;
    }
    private Vector GenerateRandomSolution(int size)
    {
      Vector v = new Vector(size);
      for (int i =0; i < v.Size; i++)
      {
        v[i] = (r.NextDouble() - 0.5) * 2;
      }
      return v;
    }

    private Matrix GenerateExpandedMatrix(Matrix source)
    {
      Matrix expMatrix = new Matrix(source.Rows + 2, source.Coloumns);
      for (int i =0; i < source.Rows; i++)
      {
        for (int j = 0; j < source.Coloumns; j++)
        {
          expMatrix[i, j] = source[i, j];
        }
      }
      for (int i =0; i < source.Coloumns; i++)
      {
        expMatrix[source.Rows, i] = source[0, i] + source[source.Rows - 1, i];
        expMatrix[source.Rows + 1, i] = source[1, i] + source[source.Rows - 2, i];
      }
      return expMatrix;
    }

    private Vector GenerateExpandedCoefs(Vector coefs)
    {
      Vector expSolution = new Vector(coefs.Size + 2);
      for (int i =0; i < coefs.Size; i++)
      {
        expSolution[i] = coefs[i];
      }
      expSolution[coefs.Size] = expSolution[0] + expSolution[coefs.Size - 1];
      expSolution[coefs.Size + 1] = expSolution[1] + expSolution[coefs.Size - 2];
      return expSolution;
    }
  }
}
