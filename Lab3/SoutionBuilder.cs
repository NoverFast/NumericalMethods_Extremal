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
      Grid = new Matrix(StepsX + 1, StepsY + 1);
      ExcatSolution = new Matrix(StepsX + 1, StepsY + 1);
      Phi = phi;
      ManagementFunc = startManagement;
      FillInitialConditions();
    }

    private void FillInitialConditions()
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
      for (int i = 0; j < Grid.Rows; i++)
      {
        Grid[i, 0] = 0;
        //ExcatSolution[i, 0] = 0;
      }
    }

    public Matrix Solve()
    {
      return null;
    }


    public Matrix SolveForward()
    {
      Array.Clear(Grid);
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
      if (writeToFile)
      {
        //Helper.ClearFile(@"G:\Универ\3 курс\Теория Разностных Схем\Лабы\2 лаба\errorImplicitValuesX" + numberOfStepsSpace + "T" + numberOfStepsTime + ".txt");
        for (int t = 0; t < gridMatrix.GetLength(0); t++)
        {
          for (int x = 0; x < gridMatrix.GetLength(1); x++)
          {
            /*  Helper.WriteToFile(@"G:\Универ\3 курс\Теория Разностных Схем\Лабы\2 лаба\errorImplicitValuesX" + numberOfStepsSpace + "T" + numberOfStepsTime + ".txt",
                t * gridStepTime, x * gridStepSpace, Math.Abs(gridMatrix[t, x] - exactSolution[t,x])); */
          }
        }
        Helper.WriteToFile(@"G:\Универ\3 курс\Теория Разностных Схем\Лабы\2 лаба\implicitUniverseMeasures.txt",
          numberOfStepsSpace, Math.Abs(Helper.UniformMeasure(gridMatrix, exactSolution)));
      }
      Console.WriteLine("Uniform Measure: " + Helper.UniformMeasure(gridMatrix, exactSolution));
      return gridMatrix;
    }

    public Matrix SolveBackwards()
    {
      return null;
    }
  }
}

 