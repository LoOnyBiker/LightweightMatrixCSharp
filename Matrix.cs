/*
    Matrix class in C#
    Written by Ivan Kuckir (ivan.kuckir@gmail.com, http://blog.ivank.net)
    Faculty of Mathematics and Physics
    Charles University in Prague
    (C) 2010
    - updated on 01.06.2014 - Trimming the string before parsing
    - updated on 14.06.2012 - parsing improved. Thanks to Andy!
    - updated on 03.10.2012 - there was a terrible bug in LU, SoLE and Inversion. Thanks to Danilo Neves Cruz for reporting that!
    - updated on 21.01.2014 - multiple changes based on comments -> see git for further info
	
    This code is distributed under MIT licence.
	
		Permission is hereby granted, free of charge, to any person
		obtaining a copy of this software and associated documentation
		files (the "Software"), to deal in the Software without
		restriction, including without limitation the rights to use,
		copy, modify, merge, publish, distribute, sublicense, and/or sell
		copies of the Software, and to permit persons to whom the
		Software is furnished to do so, subject to the following
		conditions:

		The above copyright notice and this permission notice shall be
		included in all copies or substantial portions of the Software.

		THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
		EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
		OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
		NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
		HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
		WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
		FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
		OTHER DEALINGS IN THE SOFTWARE.
*/

using System;

public partial class Matrix : BaseMatrixType
{

    Matrix L;
    Matrix U;
    Matrix P;

    private int[] permutationVector;
    private double permutationMatrixDeterminant = 1;

    public Matrix(int rows, int cols) : base(rows, cols) { }

    #region Conditions 

    public Boolean IsSquare()
    {
        return (Rows == Cols);
    }

    public Boolean IsDifferentDimentions(Matrix[] blocks)
    {
        int leftBlockRows = blocks[0].Rows + blocks[2].Rows;
        int rightBlockRows = blocks[1].Rows + blocks[3].Rows;

        int topBlockCols = blocks[0].Cols + blocks[1].Cols;
        int bottomBlockCols = blocks[2].Cols + blocks[3].Cols;

        return (leftBlockRows != rightBlockRows) || (topBlockCols != bottomBlockCols);
    }

    #endregion

    #region Base methods

    public Matrix Invert()
    {
        if (L == null) MakeLU();

        Matrix inv = new Matrix(Rows, Cols);

        for (int i = 0; i < Rows; i++)
        {
            Matrix Ei = ZeroMatrix(Rows, 1);
            Ei[i, 0] = 1;
            Matrix col = SolveWith(Ei);
            inv.SetCol(col, i);
        }
        return inv;
    }

    public double Determinant()
    {
        if (L == null) MakeLU();
        dynamic det = permutationMatrixDeterminant;
        for (int i = 0; i < Rows; i++) det *= U[i, i];
        return det;
    }

    public Matrix Transpose()
    {
        Matrix t = new Matrix(Cols, Rows);
        for (int i = 0; i < Rows; i++)
            for (int j = 0; j < Cols; j++)
                t[j, i] = this[i, j];
        return t;
    }

    #endregion

    #region Public static methods
    public static Matrix ZeroMatrix(int iRows, int iCols)
    {
        Matrix matrix = new Matrix(iRows, iCols);
        for (int i = 0; i < iRows; i++)
            for (int j = 0; j < iCols; j++)
                matrix[i, j] = (dynamic)0;
        return matrix;
    }

    public static Matrix IdentityMatrix(int iRows, int iCols)
    {
        Matrix matrix = ZeroMatrix(iRows, iCols);
        for (int i = 0; i < Math.Min(iRows, iCols); i++)
            matrix[i, i] = (dynamic)1;
        return matrix;
    }

    public static Matrix RandomMatrix(int iRows, int iCols, int dispersion)
    {
        Random random = new Random();
        Matrix matrix = new Matrix(iRows, iCols);
        for (int i = 0; i < iRows; i++)
            for (int j = 0; j < iCols; j++)
                matrix[i, j] = (dynamic)random.Next(-dispersion, dispersion);
        return matrix;
    }
    #endregion

    public Matrix[] Split(int blockRows, int blockCols)
    {
        if (blockRows > Rows || blockCols > Cols)
            throw new MException("Can't split matrix couse of wrong indexes. Checkout number of rows and cols!");

        Matrix[] result = new Matrix[4];
        Matrix A11 = new Matrix(blockRows, blockCols);
        Matrix A12 = new Matrix(blockRows, Cols - blockCols);
        Matrix A21 = new Matrix(Rows - blockRows, blockCols);
        Matrix A22 = new Matrix(Rows - blockRows, Cols - blockCols);

        #region Init A11 block
        for (int i = 0; i < blockRows; i++)
        {
            for (int j = 0; j < blockCols; j++)
            {
                A11[i, j] = this[i, j];
            }
        }
        #endregion

        #region Init A12 block
        int R12 = 0;
        for (int i = 0; i < blockRows; i++)
        {
            int C12 = 0;
            for (int j = blockCols; j < Cols; j++)
            {
                A12[R12, C12] = this[i, j];
                C12++;
            }
            R12++;
        }
        #endregion

        #region Init A21 block
        int R21 = 0;
        for (int i = blockRows; i < Rows; i++)
        {
            int C21 = 0;
            for (int j = 0; j < blockCols; j++)
            {
                A21[R21, C21] = this[i, j];
                C21++;
            }
            R21++;
        }
        #endregion

        #region Init A22 block
        int R22 = 0;
        for (int i = blockRows; i < Rows; i++)
        {
            int C22 = 0;
            for (int j = blockCols; j < Cols; j++)
            {
                A22[R22, C22] = this[i, j];
                C22++;
            }
            R22++;
        }
        #endregion

        result[0] = A11;
        result[1] = A12;
        result[2] = A21;
        result[3] = A22;
        return result;
    }
        
    public Matrix Merge(Matrix[] blocks)
    {
        int totalRows = blocks[0].Rows + blocks[2].Rows;
        int totalCols = blocks[0].Cols + blocks[3].Cols;

        if (IsDifferentDimentions(blocks))
            throw new MException("Blocks have different dimentions!");

        Matrix result = new Matrix(totalRows, totalCols);

        #region Process A11 block
        int RowLine = 0;
        int ColLine;
        for (int i = 0; i < blocks[0].Rows; i++)
        {
            ColLine = 0;
            for (int j = 0; j < blocks[0].Cols; j++)
            {
                result[RowLine, ColLine] = blocks[0][i, j];
                ColLine++;
            }
            RowLine++;
        }
        #endregion

        #region Process A12 block
        RowLine = 0;
        ColLine = blocks[0].Cols;
        for (int i = 0; i < blocks[1].Rows; i++)
        {
            ColLine = blocks[0].Cols;
            for (int j = 0; j < blocks[1].Cols; j++)
            {
                result[RowLine, ColLine] = blocks[1][i, j];
                ColLine++;
            }
            RowLine++;
        }
        #endregion

        #region Process A21 block
        RowLine = blocks[0].Rows;
        ColLine = 0;
        for (int i = 0; i < blocks[2].Rows; i++)
        {
            ColLine = 0;
            for (int j = 0; j < blocks[2].Cols; j++)
            {
                result[RowLine, ColLine] = blocks[2][i, j];
                ColLine++;
            }
            RowLine++;
        }
        #endregion

        #region Process A22 block
        RowLine = blocks[0].Rows;
        ColLine = blocks[0].Cols;
        for (int i = 0; i < blocks[3].Rows; i++)
        {
            ColLine = blocks[0].Cols;
            for (int j = 0; j < blocks[3].Cols; j++)
            {
                result[RowLine, ColLine] = blocks[3][i, j];
                ColLine++;
            }
            RowLine++;
        }
        #endregion

        return result;
    }

    public Matrix GetCol(int k)
    {
        k--;
        if (isColExist(k))
        {
            Matrix m = new Matrix(Rows, 1);
            for (int i = 0; i < Rows; i++)
                m[i, 0] = this[i, k];
            return m;
        }
        else
            throw new MException("Dimention");
    }

    public void SetCol(Matrix vector, int k)
    {
        k--;
        if (isColExist(k))
        {            
            for (int i = 0; i < Rows; i++)
                this[i, k] = vector[i, 0];
        }
        else
            throw new MException("Dimention!");        
    }

    public void RowAddition(int toIndex, int fromIndex)
    {
        for (int i = 0; i < Cols; i++)
        {
            this[toIndex, i] += this[fromIndex, i];
        }
    }

    public void RowMultiplication(int index, double value)
    {
        for (int i = 0; i < Cols; i++)
        {
            this[index, i] *= value;
        }
    }

    public void RowMultiplicationAndAddition(int toIndex, int fromIndex, double value)
    {
        RowMultiplication(fromIndex, value);
        RowAddition(toIndex, fromIndex);
    }

    public void RowSwap(int firstIndex, int secondIndex)
    {
        for (int i = 0; i < Cols; i++)
        {
            double tmp = this[firstIndex, i];
            this[firstIndex, i] = this[secondIndex, i];
            this[secondIndex, i] = tmp;
        }
    }

    public void ColSwap(int firstIndex, int secondIndex)
    {
        for (int i = 0; i < Cols; i++)
        {
            double tmp = this[i, firstIndex];
            this[i, firstIndex] = this[secondIndex, i];
            this[i, secondIndex] = tmp;
        }
    }

    public void ColDivide(int index, double value)
    {
        for (int i = 0; i < Cols; i++)
        {
            this[i, index] /= value;
        }
    }

    public Matrix GetPermutationMatrix()
    {
        if (L == null) MakeLU();

        Matrix matrix = ZeroMatrix(Rows, Cols);
        for (int i = 0; i < Rows; i++)
            matrix[permutationVector[i], i] = (dynamic)1;
        return matrix;
    }

    public Matrix Power(int pow)           // Power matrix to exponent
    {
        if (pow == 0) return IdentityMatrix(Rows, Cols);
        if (pow == 1) return Duplicate();
        if (pow == -1) return Invert();

        Matrix x;
        if (pow < 0) { x = Invert(); pow *= -1; }
        else x = Duplicate();

        Matrix ret = IdentityMatrix(Rows, Cols);
        while (pow != 0)
        {
            if ((pow & 1) == 1) ret *= x;
            x *= x;
            pow >>= 1;
        }
        return ret;
    }

    public Matrix Duplicate()
    {
        Matrix matrix = new Matrix(Rows, Cols);
        for (int i = 0; i < Rows; i++)
            for (int j = 0; j < Cols; j++)
                matrix[i, j] = this[i, j];
        return matrix;
    }
}

//  The class for exceptions

public class MException : Exception
{
    public MException(string Message)
        : base(Message)
    { }
}