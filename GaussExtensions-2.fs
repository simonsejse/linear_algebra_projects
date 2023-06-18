module ProjectB
// LinalgDat23
// Authors: Francois Lauze

open System
open LinAlgDat.Core

type GaussOps =
    class

        /// <summary>
        /// This function creates an augmented Matrix given a Matrix A, and a
        /// right-hand side Vector v.
        /// </summary>
        ///
        /// <remarks>
        /// See page 12 in "Linear Algebra for Engineers and Scientists"
        /// by K. Hardy.
        /// </remarks>
        ///
        /// <param name="A">An M-by-N Matrix.</param>
        /// <param name="v">An M-size Vector.</param>
        ///
        /// <returns>An M-by-(N+1) augmented Matrix [A | v].</returns>
        static member AugmentRight (A: Matrix) (v: Vector) : Matrix =
            let m_rows: int = A.M_Rows
            let n_cols: int = A.N_Cols

            let retval: float array2d = Array2D.zeroCreate m_rows (n_cols + 1)

            for i: int32 in 0 .. m_rows - 1 do
                for j: int32 in 0 .. n_cols - 1 do
                    retval.[i, j] <- A.[i, j]

                retval.[i, n_cols] <- v.[i]

            Matrix retval

        /// <summary>
        /// This function computes the elementary row replacement operation on
        /// the given matrix.
        /// </summary>
        ///
        /// <remarks>
        /// Note that we add the row (as in the lectures) instead of subtracting
        /// the row (as in the textbook).
        /// </remarks>
        ///
        /// <param name="A">
        /// An M-by-N matrix to perform the elementary row operation on.
        /// </param>
        /// <param name="i">
        /// The index of the row to replace.
        /// </param>
        /// <param name="m">
        /// The multiple of row j to add to row i.
        /// </param>
        /// <param name="j">
        /// The index of the row whose mutiple is added to row i.
        /// </param>
        ///
        /// <returns>
        /// The resulting M-by-N matrix after having performed the elementary
        /// row operation.
        /// </returns>
        static member ElementaryRowReplacement (A: Matrix) (i: int) (m: float) (j: int) : Matrix =
            //The elementary row operation: mj + i -> i
            let retval: Matrix = Matrix A

            for column in 0 .. retval.N_Cols - 1 do
                retval.[i, column] <- (retval.Row(j) * m).[column] + retval.Row(i).[column]

            retval

        /// <summary>
        /// This function computes the elementary row interchange operation on
        /// the given matrix.
        /// </summary>
        ///
        /// <param name="A">
        /// An M-by-N matrix to perform the elementary row operation on.
        /// </param>
        /// <param name="i">
        /// The index of the first row of the rows to interchange.
        /// </param>
        /// <param name="j">
        /// The index of the second row of the rows to interchange.
        /// </param>
        ///
        /// <returns>
        /// The resulting M-by-N matrix after having performed the elementary
        /// row operation.
        /// </returns>
        static member ElementaryRowInterchange (A: Matrix) (i: int) (j: int) : Matrix =
            let retval: Matrix = Matrix A

            let tempIRow: Vector = retval.Row(i)

            for column: int32 in 0 .. retval.N_Cols - 1 do
                retval.[i, column] <- retval.Row(j).[column]
                retval.[j, column] <- tempIRow.[column]

            retval

        /// <summary>
        /// This function computes the elementary row scaling operation on the
        /// given matrix.
        /// </summary>
        ///
        /// <param name="A">
        /// An M-by-N matrix to perform the elementary row operation on.
        /// </param>
        /// <param name="i">The index of the row to scale.</param>
        /// <param name="c">The value to scale the row by.</param>
        ///
        /// <returns>
        /// The resulting M-by-N matrix after having performed the elementary
        /// row operation.
        /// </returns>
        static member ElementaryRowScaling (A: Matrix) (i: int) (c: float) : Matrix =
            let retval: Matrix = A

            for column: int32 in 0 .. retval.N_Cols - 1 do
                retval[i, column] <- retval.Row(i).[column] * c

            retval

        /// <summary>
        /// This function executes the forward reduction algorithm provided in
        /// the assignment text to achieve row Echelon form of a given
        /// augmented matrix.
        /// </summary>
        ///
        /// <param name="A">
        /// An M-by-N matrix, augmented (or not).
        /// </param>
        ///
        /// <returns>
        /// An M-by-N matrix that is the row Echelon form.
        /// </returns>

        static member ForwardReduction(A: Matrix) : Matrix =
            let mutable retval = new Matrix(A)

            let tolerance = 0.00000001

            let M = A.M_Rows
            let N = A.N_Cols

            let mutable i = 0
            let mutable j = 0

            while i < M && j < N do
                if abs (retval.[i, j]) < tolerance then
                    let mutable k = i + 1

                    while k < M && abs (retval.[k, j]) < tolerance do
                        k <- k + 1

                    if k < M then
                        retval <- GaussOps.ElementaryRowInterchange retval i k

                if abs (retval.[i, j]) > tolerance then
                    let mutable k = i + 1

                    while k < M do
                        retval <- GaussOps.ElementaryRowReplacement retval k (-retval.[k, j] / retval.[i, j]) i
                        k <- k + 1

                    i <- i + 1

                j <- j + 1

            retval

        /// <summary>
        /// This function executes the backward reduction algorithm provided in
        /// the assignment text given an augmented matrix in row Echelon form.
        /// </summary>
        ///
        /// <param name="A">
        /// An M-by-N augmented matrix in row Echelon form.
        /// </param>
        ///
        /// <returns>
        /// The resulting M-by-N matrix after executing the algorithm.
        /// </returns>
        static member BackwardReduction(A: Matrix) : Matrix =
            let tolerance = 0.00000001

            let mutable U = new Matrix(A)

            let M = A.M_Rows
            let N = A.N_Cols

            let mutable i = M - 1
            let mutable j = N - 1

            while i >= 0 && j >= 0 do
                let mutable k = M - 1

                while k >= 0 do
                    if abs (U.[k, j]) > tolerance then
                        U <- GaussOps.ElementaryRowScaling U k (1.0 / U.[k, j])

                        let mutable l = k - 1

                        while l >= 0 do
                            U <- GaussOps.ElementaryRowReplacement U l (-U.[l, j]) k
                            l <- l - 1

                    i <- i + 1
                    k <- k - 1

                j <- j - 1
                i <- i - 1

            U


        /// <summary>
        /// This function performs Gauss elimination of a linear system
        /// given in matrix form by a coefficient matrix and a right hand side
        /// vector. It is assumed that the corresponding linear system is
        /// consistent and has exactly one solution.
        /// </summary>
        ///
        /// <remarks>
        /// Hint: Combine ForwardReduction and BackwardReduction.
        /// </remarks>
        ///
        /// <param name="A">An M-by-N Matrix.</param>
        /// <param name="b">An M-size Vector.</param>
        ///
        /// <returns>The N-sized vector x such that A * x = b.</returns>
        static member GaussElimination (A: Matrix) (b: Vector) : Vector =
            let ExtractX (matrix: Matrix) : Vector =
                let mutable retval = new Vector(matrix.M_Rows)

                for row in 0 .. matrix.M_Rows - 1 do
                    retval.[row] <- matrix.[row, matrix.N_Cols - 1]

                retval

            b
            |> GaussOps.AugmentRight A
            |> GaussOps.ForwardReduction
            |> GaussOps.BackwardReduction
            |> ExtractX





































































































    end
