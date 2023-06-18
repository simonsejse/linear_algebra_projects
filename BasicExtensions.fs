module ProjectA

open System
open LinAlgDat.Core

type BasicOps =
    class
        /// <summary>
        /// This function creates an augmented Matrix given a Matrix A, and a
        /// right-hand side Vector v.
        /// </summary>
        ///
        /// <remarks>
        /// See page 12 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// This implementation is provided for you.
        /// </remarks>
        ///
        /// <param name="A">An M-by-N Matrix.</param>
        /// <param name="v">An M-size Vector.</param>
        ///
        /// <returns>An M-by-(N+1) augmented Matrix [A | v].</returns>
        static member AugmentRight (A: Matrix) (v: Vector) : Matrix =
            let m_rows = A.M_Rows
            let n_cols = A.N_Cols

            let retval = Array2D.zeroCreate m_rows (n_cols + 1)

            for i in 0 .. m_rows - 1 do
                for j in 0 .. n_cols - 1 do
                    retval.[i, j] <- A.[i, j]

                retval.[i, n_cols] <- v.[i]

            Matrix retval

        /// <summary>
        /// This function computes the Matrix-Vector product of a Matrix A,
        /// and a column Vector v.
        /// </summary>
        ///
        /// <remarks>
        /// See page 68 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="A">An M-by-N Matrix.</param>
        /// <param name="v">An N-size Vector.</param>
        ///
        /// <returns>An M-size Vector b such that b = A * v.</returns>
        static member MatVecProduct (A: Matrix) (v: Vector) : Vector =
            let mSize = A.M_Rows //An M-size Vector b
            let mutable b = Vector(mSize)

            //First we will iterate through each row
            for row in 0 .. A.M_Rows - 1 do
                //Then we will iterate through each column
                for column in 0 .. A.N_Cols - 1 do
                    //Nu sætter vi b = A*v
                    b.[row] <- b.[row] + A.[row, column] * v.[column]

            b

        /// <summary>
        /// This function computes the Matrix product of two given matrices
        /// A and B.
        /// </summary>
        ///
        /// <remarks>
        /// See page 58 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarls>
        ///
        /// <param name="A">An M-by-N Matrix.</param>
        /// <param name="B">An N-by-P Matrix.</param>
        ///
        /// <returns>The M-by-P Matrix A * B.</returns>
        static member MatrixProduct (A: Matrix) (B: Matrix) : Matrix =
            let m = A.M_Rows //Rows in A
            let p = B.N_Cols //Columns in B
            let mutable C = Matrix(m, p)
            //Iterate through each row in A
            for row in 0 .. m - 1 do
                //Iterate through each column in B
                for column in 0 .. p - 1 do
                    //For each column value we want to multiply
                    // the row values in A with the corresponding column values in B
                    //We can do this with another for loop k that goes from 0
                    //to the length of the column, and then just multiply the values
                    //together and add them to the sum and store in the new Matrix C.
                    for k in 0 .. A.N_Cols - 1 do
                        C.[row, column] <- C.[row, column] + A.[row, k] * B.[k, column]

            C

        /// <summary>
        /// This function computes the transpose of a given Matrix.
        /// </summary>
        ///
        /// <remarks>
        /// See page 69 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="A">An M-by-N Matrix.</param>
        ///
        /// <returns>The N-by-M Matrix B such that B = A^T.</returns>
        static member Transpose(A: Matrix) : Matrix =
            let M = A.M_Rows
            let N = A.N_Cols
            let mutable B = Matrix(N, M)

            //Iterate through each row in A
            for row in 0 .. M - 1 do
                //Just swap column and row around
                for column in 0 .. N - 1 do
                    B.[column, row] <- A.[row, column]

            B

        /// <summary>
        /// This function computes the Euclidean Vector norm of a given
        /// Vector.
        /// </summary>
        ///
        /// <remarks>
        /// See page 197 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="v">An N-dimensional Vector.</param>
        ///
        /// <returns>The Euclidean norm of the Vector.</returns>
        static member VectorNorm(v: Vector) : float =
            let mutable sum: float = 0.0
            //Pretty straight forward:
            //In the book page 197 it states the formula
            //is the sum of all vector coordinates squared
            //then the whole sum square rooted:
            for coord in 0 .. v.Size - 1 do
                sum <- sum + v.[coord] ** 2.0

            sum |> sqrt









































































    end
