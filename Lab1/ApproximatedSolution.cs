using System;
using MathPrimitiveLibrary;

namespace AFLETUNOV_LR1.Lab1
{
  public class ApproximatedSolution
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
    public Matrix ApproximatedCoefs { get; private set; }

    public Random r { get; private set; }

    // в качестве начального приближения метода ньютона для решения СНАУ можно использовать точное решение
    public ApproximatedSolution(int size)
    {
      InitializeData(size); // Инициализируем все начальные данные, с которыми впоследствие будем работать
    }

    private void InitializeData(int size)
    {
      SourceMatrix = GenerateRandomMatrix(size);
      ExactSolution = GenerateRandomSolution(size);
      FreeCoefs = SourceMatrix * ExactSolution;
      ExpandedMatrix = GenerateExpandedMatrix(SourceMatrix);
      ExpandedCoefs = GenerateExpandedCoefs(FreeCoefs);
      ApproximatedMatrix = AddNoise(ExpandedMatrix);
      //test
    }

    private Matrix AddNoise(Matrix m)
    {
      Random r = new Random();
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
    private Vector AddNoise(Vector v)
    {
      Random r = new Random();
      Vector noiseVector = new Vector(v);
      for (int i = 0; i < noiseVector.Size; i++)
      {
        noiseVector[i] *= (1 + ((r.NextDouble() - 0.5) * 2) * 0.01);
      }
      return noiseVector;
    }
    private Matrix GenerateRandomMatrix(int size)
    {
      Random r = new Random();
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
      Random r = new Random();
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
      for (int i =0; i < source.Coloumns; i++)
      {
        source[10, i] = source[0, i] + source[9, i];
        source[11, i] = source[1, i] + source[8, i];
      }
      return expMatrix;
    }

    private Vector GenerateExpandedCoefs(Vector coefs)
    {
      Vector expSolution = new Vector(coefs.Size + 2);
      expSolution[10] = expSolution[0] + expSolution[9];
      expSolution[11] = expSolution[1] + expSolution[8];
      return expSolution;
    }
  }
}
