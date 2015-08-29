using System;
using System.Text;

public class BaseMatrixType : ICloneable
{

    public int Rows { get; protected set; }
    public int Cols { get; protected set; }
    protected double[,] items;

    public double this[int row, int col]
    {
        get { return items[row, col]; }
        set { items[row, col] = value; }
    }

    public BaseMatrixType(int rows, int cols)
    {
        Rows = rows;
        Cols = cols;
        items = new double[rows, cols];
    }

    #region Conditions
    public bool isRowExist(int index)
    {
        return index >= 0 && index < Rows;
    }

    public bool isColExist(int index)
    {
        return index >= 0 && index < Cols;
    }

    public bool IsMultiplicationPossible(Matrix other)
    {
        return Cols == other.Rows;
    }
    #endregion

    protected static int GetMaxDimention(BaseMatrixType m1, BaseMatrixType m2)
    {
        int first = Math.Max(m1.Rows, m1.Cols);
        int second = Math.Max(m2.Rows, m2.Cols);
        return Math.Max(first, second);
    }

    protected int GetMaxDimention(BaseMatrixType m)
    {
        return Math.Max(m.Rows, m.Cols);
    }

    #region Operators
    protected static BaseMatrixType Zero(int iRows, int iCols)
    {
        BaseMatrixType matrix = new BaseMatrixType(iRows, iCols);
        for (int i = 0; i < iRows; i++)
            for (int j = 0; j < iCols; j++)
                matrix[i, j] = 0;
        return matrix;
    }

    private static BaseMatrixType Multiply(BaseMatrixType m1, BaseMatrixType m2)
    {
        if (m1.Cols != m2.Rows) throw new MException("Wrong dimensions of matrix!");

        BaseMatrixType result = Zero(m1.Rows, m2.Cols);
        for (int i = 0; i < result.Rows; i++)
            for (int j = 0; j < result.Cols; j++)
                for (int k = 0; k < m1.Cols; k++)
                    result[i, j] += m1[i, k] * m2[k, j];
        return result;
    }

    private static BaseMatrixType Multiply(double n, BaseMatrixType m)
    {
        BaseMatrixType r = new BaseMatrixType(m.Rows, m.Cols);
        for (int i = 0; i < m.Rows; i++)
            for (int j = 0; j < m.Cols; j++)
                r[i, j] = m[i, j] * n;
        return r;
    }

    private static BaseMatrixType Add(BaseMatrixType m1, BaseMatrixType m2)         // Sčítání matic
    {
        if (m1.Rows != m2.Rows || m1.Cols != m2.Cols) throw new MException("Matrices must have the same dimensions!");
        BaseMatrixType r = new BaseMatrixType(m1.Rows, m1.Cols);
        for (int i = 0; i < r.Rows; i++)
            for (int j = 0; j < r.Cols; j++)
                r[i, j] = m1[i, j] + m2[i, j];
        return r;
    }

    public static BaseMatrixType operator -(BaseMatrixType m)
    { return Multiply(-1, m); }

    public static BaseMatrixType operator +(BaseMatrixType m1, BaseMatrixType m2)
    { return Add(m1, m2); }

    public static BaseMatrixType operator -(BaseMatrixType m1, BaseMatrixType m2)
    { return Add(m1, -m2); }

    public static BaseMatrixType operator *(BaseMatrixType m1, BaseMatrixType m2)
    { return Multiply(m1, m2); }

    public static BaseMatrixType operator *(double n, BaseMatrixType m)
    { return Multiply(n, m); }
    #endregion

    #region IClonable
    private BaseMatrixType Duplicate()
    {
        BaseMatrixType matrix = new BaseMatrixType(Rows, Cols);
        for (int i = 0; i < Rows; i++)
            for (int j = 0; j < Cols; j++)
                matrix[i, j] = this[i, j];
        return matrix;
    }

    public object Clone()
    {
        return Duplicate();
    }
    #endregion

    public override string ToString()                           // Function returns matrix as a string
    {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i < Rows; i++)
        {
            for (int j = 0; j < Cols; j++)
                s.Append(String.Format("{0:0.00}", this[i, j]) + "\t");
            s.AppendLine();
        }
        return s.ToString();
    }
}
