using System;
using MathPrimitivesLibrary;
using System.IO;

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
    public Vector ExpandedSolution { get; private set; }

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

      /*SourceMatrix = new Matrix(new double[,] {
        { 0.700976, 0.809676, 0.0887955, 0.121479},
        { 0.348307, 0.421962, 0.699805, 0.0663843},
        { 0.587482, 0.642966, 0.990603, 0.295718},
        { 0.271337, 0.069656, 0.949639, 0.382175}
      });
      ExactSolution = new Vector(new double[] { 0.32271, 0.988366, 0.842522, 0.275314 });
      FreeCoefs = new Vector(new double[] { 1.13473, 1.13733, 1.74109, 1.06172 });
      //FreeCoefs = SourceMatrix * ExactSolution;

      ExpandedMatrix = new Matrix(new double[,] {
        { 0.700976, 0.809676, 0.0887955, 0.121479},
        { 0.348307, 0.421962, 0.699805, 0.0663843},
        { 0.587482, 0.642966, 0.990603, 0.295718},
        { 0.271337, 0.069656, 0.949639, 0.382175},
        { 0.972313, 0.879332, 1.03843, 0.503654 },
        { 0.935789, 1.06493, 1.69041, 0.362103 }
      });
      ExpandedSolution = new Vector(new double[] { 0.32271, 0.988366, 0.842522, 0.275314, 0.598024, 1.83089 });
      ExpandedCoefs = new Vector(new double[] { 1.13473, 1.13733, 1.74109, 1.06172, 2.19644, 2.87843 });

      ApproximatedMatrix = new Matrix(new double[,] {
        { 0.70099, 0.809685, 0.0888016, 0.121472},
        { 0.348334, 0.421992, 0.69978, 0.0663791},
        { 0.587499, 0.643005, 0.990674, 0.295732},
        { 0.271321, 0.0696512, 0.949663, 0.382158},
        { 0.972268, 0.87929, 1.03834, 0.503658 },
        { 0.935881, 1.06502, 1.69046, 0.362137 }
      });

      ApproximatedCoefs = new Vector(new double[] { 1.1347, 1.13723, 1.74124, 1.06177, 2.19637, 2.87858 });

      Console.WriteLine(); */
      SourceMatrix = GenerateRandomMatrix(size);
      ExactSolution = GenerateRandomSolution(size);
      FreeCoefs = SourceMatrix * ExactSolution;
      ExpandedMatrix = GenerateExpandedMatrix(SourceMatrix);
      ExpandedCoefs = GenerateExpandedCoefs(FreeCoefs);
      ApproximatedMatrix = AddNoise(ExpandedMatrix, 1e-5);
      ApproximatedCoefs = AddNoise(ExpandedCoefs, 1e-5);

#if DEBUG
      Console.WriteLine("\n\tSource Matrix | Coeffs");
      SourceMatrix.Show();
      FreeCoefs.Show();

      Console.WriteLine("\n\tExpanded Matrix | Expanded Coeffs");
      ExpandedMatrix.Show();
      ExpandedCoefs.Show();

      Console.WriteLine("\n\tExpanded Matrix (With Noise) | Expanded Coeffs (Without Noise)");
      ApproximatedMatrix.Show();
      ApproximatedCoefs.Show(); 
#endif
    }

    private void ClearFile(string path)
    {
      File.WriteAllText(path, string.Empty);
    }

    public void Calculate()
    {
      SigmaAndH();
      FindLambdaDelta();
      FindSolution();
    }
    double LFunc(Vector u)
    {
      double buff0, buff1 = 0, buff2 = 0;
      for (int i =0; i <  ExpandedCoefs.Size; i++)
      {
        buff0 = 0;
        for (int j = 0; j < FreeCoefs.Size; j++)
        {
          buff0 += ApproximatedMatrix[i, j] * u[j]; 
        }
        buff0 -= ApproximatedCoefs[i];
        buff1 += Math.Pow(buff0, 2);
      }
      for (int i =0; i < FreeCoefs.Size; i++)
      {
        buff2 += Math.Pow(u[i], 2);
      }
      //Console.WriteLine("---");
      //Console.WriteLine(Math.Sqrt(buff1) + H * Math.Sqrt(buff2) + Sigma);
      //Console.WriteLine("---");
      return Math.Sqrt(buff1) + H * Math.Sqrt(buff2) + Sigma;
    }

    private double dL(Vector u, int k)
    {
      double buff0, buff1 = 0, buff2 = 0, buff3 = 0;
      for (int i = 0; i < ExpandedCoefs.Size; i++)
      {
        buff0 = 0;
        for (int j =0; j < FreeCoefs.Size; j++)
        {
          buff0 += ApproximatedMatrix[i, j] * u[j];
        }
        buff0 -= ApproximatedCoefs[i];
        buff1 += ApproximatedMatrix[i, k] * buff0;
        buff2 += Math.Pow(buff0, 2);
      }
      for (int i =0; i < FreeCoefs.Size; i++)
      {
        buff3 += Math.Pow(u[i], 2);
      }
      return buff1 / Math.Sqrt(buff2) + H * u[k] / Math.Sqrt(buff3);
    }

    private Vector gradL(Vector u)
    {
      Vector gradient = new Vector(ExactSolution.Size);
      for (int i =0; i < gradient.Size; i++)
      {
        gradient[i] = dL(u, i);
        //Console.WriteLine(gradient[i]);
      }
      return gradient; 
    }

    private Vector CGM(Matrix A, Vector b, double precision)
    {
      int n = A.Rows;
      Vector xPrev = new Vector(n);
      Vector x = new Vector(xPrev.Size);
      Vector rPrev = b - A * xPrev;

      Vector r = new Vector(rPrev.Size);
      rPrev.CopyTo(r);
      Vector zPrev = new Vector(rPrev.Size);
      rPrev.CopyTo(zPrev);
      Vector z = new Vector(zPrev.Size);
      double error;
      do
      {
        double alphaK = rPrev.DotProduct(rPrev) / (A * zPrev).DotProduct(zPrev);
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
        for (int j =0; j < SourceMatrix.Rows; j++)
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
      Console.WriteLine("\t--- Sigma | H ---");
      Matrix aTmp = ExpandedMatrix - ApproximatedMatrix;
      Vector fTmp = ExpandedCoefs - ApproximatedCoefs;
      H = MatrixNorm(aTmp);
      Sigma = fTmp.Norm();
      Console.WriteLine($"H: {H}, | Sigma: {Sigma}");

    }

    private double MatrixNorm(Matrix a)
    {
      double sum = 0;
      for (int i =0; i < a.Rows;i++)
      {
        for (int j =0; j < a.Coloumns;j++)
        {
          sum += Math.Pow(a[i, j], 2);
        }
      }
      return Math.Sqrt(sum);
    }

    private Vector GradientDescent(Func<Vector, double> f, Func<Vector, Vector> gradF, Vector x0, double precision, int M )
    {
      int k = 0;
      Vector xPrev = new Vector(x0);
      Vector x = new Vector(xPrev.Size);
      double tK = 3;
      //gradF(xPrev).Show();
      while (true)
      {
        k += 1;
        Vector gradientX = new Vector(xPrev.Size);
        gradientX = gradF(xPrev);
        for (int i =0; i < gradientX.Size; i++)
        {
          gradientX[i] *= tK;
        }
        (xPrev - gradientX).CopyTo(x);
        /*Console.WriteLine("\nxPrev:");
        xPrev.Show();
        Console.WriteLine("\ngradient:");
        gradientX.Show();
        Console.WriteLine("\nx:");
        x.Show();
        Console.WriteLine("\nf(x): {0}", f(x));
        Console.WriteLine("f(xPrev): {0}", f(xPrev));*/
        if (f(x) - f(xPrev) < 0)
        {
          Vector tmp = x - xPrev;
          //tmp.Show();
          if (tmp.Norm() <= precision || k + 1 > M)
          {
            return x;
          }
        }
        else
        {
          tK /= 2;
        }
        for (int i =0; i < xPrev.Size; i++)
        {
          xPrev[i] = x[i];
        }
      }
    }

    private void FindLambdaDelta()
    {
      Console.WriteLine("\t--- Lambda-Delta ---");
      Vector x0 = new Vector(ExactSolution);
      Console.WriteLine("x0: ");
      x0.Show();
      Vector LMin = GradientDescent(LFunc, gradL, x0, 1e-6, 1000);
      LambdaDelta = LFunc(LMin);
      Console.WriteLine("\nLMin: ");
      LMin.Show();
      Console.WriteLine($"\nLambda Delta: {LambdaDelta}");
    }
    private void FindSolution()
    {
      Console.WriteLine("\t--- Solution ---");
      double alpha = 1;
      Matrix matrix = ApproximatedMatrix.Transpose() * ApproximatedMatrix + SourceMatrix.IdentityMatrix * alpha;
      Vector vector = ApproximatedMatrix.Transpose() * ApproximatedCoefs;
      Vector answer = CGM(matrix, vector, 1e-6);

      ClearFile("../../Lab1/Results/rho.txt");
      StreamWriter rho = new StreamWriter("../../Lab1/Results/rho.txt");
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
        //Console.WriteLine(alpha + " " + val);
        if (alpha <= 0.25)
        {
          rho.WriteLine(alpha + " " + val);
        }
        if (val < 1e-4 && alpha >= 0)
        {
          needAlpha = alpha;
          needValue = val;
          needAnswer = answer;
        }
      }
      rho.Close();
      Console.WriteLine($"\nValue: {needValue}");
      Console.WriteLine($"Alpha: {needAlpha}");
      Console.WriteLine("Approximated Solution: ");
      needAnswer.Show();
      Console.WriteLine($"\nAnswer Norm: {needAnswer.Norm()}");
      Console.WriteLine("\nExact Solution: ");
      ExactSolution.Show();
      Console.WriteLine($"\nExact Solution Norm: {ExactSolution.Norm()}");
      Console.WriteLine($"Norm Difference: {needAnswer.Norm() - ExactSolution.Norm()}");
    }

    private Matrix AddNoise(Matrix m, double prec)
    {
      Matrix noiseMatrix = new Matrix(m);
      for (int i = 0; i < noiseMatrix.Rows; i++)
      {
        for (int j = 0; j < noiseMatrix.Coloumns; j++)
        {
          noiseMatrix[i, j] *= (1 + ((r.NextDouble() - 0.5) * 2) * prec);
        }
      }
      return noiseMatrix;
    }
    private Vector AddNoise(Vector v, double prec)
    {
      Vector noiseVector = new Vector(v);
      for (int i = 0; i < noiseVector.Size; i++)
      {
        noiseVector[i] *= (1 + ((r.NextDouble() - 0.5) * 2) * prec);
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
          randomValue = Math.Abs((r.NextDouble() - 0.5) * 2);
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
        v[i] = Math.Abs((r.NextDouble() - 0.5) * 2);
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
