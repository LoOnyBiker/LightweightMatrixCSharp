using System;


public partial class Matrix
{

    public void MakeLU()                        // Function for LU decomposition
    {
        if (!IsSquare()) throw new MException("The matrix is not square!");
        L = IdentityMatrix(Rows, Cols);
        U = Duplicate();

        permutationVector = new int[Rows];
        for (int i = 0; i < Rows; i++) permutationVector[i] = i;

        double p = 0;
        double pom2;
        int k0 = 0;
        int pom1 = 0;

        for (int k = 0; k < Cols - 1; k++)
        {
            p = 0;
            for (int i = k; i < Rows; i++)      // find the row with the biggest pivot
            {
                if (Math.Abs(U[i, k]) > p)
                {
                    p = Math.Abs(U[i, k]);
                    k0 = i;
                }
            }
            if (p == 0) // samé nuly ve sloupci
                throw new MException("The matrix is singular!");

            pom1 = permutationVector[k];
            permutationVector[k] = permutationVector[k0];   // switch two rows in permutation matrix
            permutationVector[k0] = pom1;

            for (int i = 0; i < k; i++)
            {
                pom2 = L[k, i]; L[k, i] = L[k0, i]; L[k0, i] = pom2;
            }

            if (k != k0) permutationMatrixDeterminant *= -1;

            for (int i = 0; i < Cols; i++)                  // Switch rows in U
            {
                pom2 = U[k, i]; U[k, i] = U[k0, i]; U[k0, i] = pom2;
            }

            for (int i = k + 1; i < Rows; i++)
            {
                L[i, k] = U[i, k] / U[k, k];
                for (int j = k; j < Cols; j++)
                    U[i, j] = U[i, j] - L[i, k] * U[k, j];
            }
        }
    }

    public static Matrix SubsForth(Matrix A, Matrix b)          // Function solves Ax = b for A as a lower triangular matrix
    {
        if (A.L == null) A.MakeLU();
        int n = A.Rows;
        Matrix x = new Matrix(n, 1);

        for (int i = 0; i < n; i++)
        {
            x[i, 0] = b[i, 0];
            for (int j = 0; j < i; j++) x[i, 0] -= A[i, j] * x[j, 0];
            x[i, 0] = x[i, 0] / A[i, i];
        }
        return x;
    }

    public static Matrix SubsBack(Matrix A, Matrix b)           // Function solves Ax = b for A as an upper triangular matrix
    {
        if (A.L == null) A.MakeLU();
        int n = A.Rows;
        Matrix x = new Matrix(n, 1);

        for (int i = n - 1; i > -1; i--)
        {
            x[i, 0] = b[i, 0];
            for (int j = n - 1; j > i; j--) x[i, 0] -= A[i, j] * x[j, 0];
            x[i, 0] = x[i, 0] / A[i, i];
        }
        return x;
    }

    public Matrix SolveWith(Matrix v)                        // Function solves Ax = v in confirmity with solution vector "v"
    {
        if (Rows != Cols) throw new MException("The matrix is not square!");
        if (Rows != v.Rows) throw new MException("Wrong number of results in solution vector!");
        if (L == null) MakeLU();

        Matrix b = new Matrix(Rows, 1);
        for (int i = 0; i < Rows; i++) b[i, 0] = v[permutationVector[i], 0];   // switch two items in "v" due to permutation matrix

        Matrix z = SubsForth(L, b);
        Matrix x = SubsBack(U, z);

        return x;
    }

    // TODO check for redundancy with MakeLU() and SolveWith()
    public void MakeRref()                                    // Function makes reduced echolon form
    {
        int lead = 0;
        for (int r = 0; r < Rows; r++)
        {
            if (Cols <= lead) break;
            int i = r;
            while (this[i, lead] == 0)
            {
                i++;
                if (i == Rows)
                {
                    i = r;
                    lead++;
                    if (Cols == lead)
                    {
                        lead--;
                        break;
                    }
                }
            }
            for (int j = 0; j < Cols; j++)
            {
                double temp = this[r, j];
                this[r, j] = this[i, j];
                this[i, j] = temp;
            }
            double div = this[r, lead];
            for (int j = 0; j < Cols; j++) this[r, j] /= div;
            for (int j = 0; j < Rows; j++)
            {
                if (j != r)
                {
                    double sub = this[j, lead];
                    for (int k = 0; k < Cols; k++) this[j, k] -= (sub * this[r, k]);
                }
            }
            lead++;
        }
    }

}
