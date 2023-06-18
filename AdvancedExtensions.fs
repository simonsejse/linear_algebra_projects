module ProjectC

open System
open LinAlgDat.Core

type AdvancedOps = class

    /// <summary>
    ///     This function creates the square submatrix given a square matrix as
    ///     well as row and column indices to remove from it.
    /// </summary>
    /// <remarks>
    ///     See page 246-247 in "Linear Algebra for Engineers and Scientists"
    ///     by K. Hardy.
    /// </remarks>
    /// <param name="A">An N-by-N matrix.</param>
    /// <param name="i">The index of the row to remove.</param>
    /// <param name="j">The index of the column to remove.</param>
    /// <returns>The resulting (N - 1)-by-(N - 1) submatrix.</returns>
    static member SquareSubMatrix (A : Matrix) (i : int) (j : int) : Matrix =
       let N = A.N_Cols - 1

       let mutable M = Matrix(N, N)
       
       let mutable rowOffset = 0
       let mutable colOffset = 0

       for row = 0 to N - 1 do
           if row = i then
               rowOffset <- 1
               colOffset <- 0

           for col = 0 to N - 1 do
               if col = j then
                   colOffset <- 1
               M[row, col] <- A[row + rowOffset, col + colOffset]
           colOffset <- 0
       M    
              
    /// <summary>
    ///     This function computes the determinant of a given square matrix.
    /// </summary>
    /// <remarks>
    ///     See page 247 in "Linear Algebra for Engineers and Scientists"
    ///     by K. Hardy.
    /// </remarks>
    /// <remarks>
    ///     Hint: Use SquareSubMatrix.
    /// </remarks>
    /// <param name="A">An N-by-N matrix.</param>
    /// <returns>The determinant of the matrix</returns>
    static member Determinant (A : Matrix) : float =
        //Using the recursive definition for determinants
        let rec det (subA:Matrix) =
                if (subA.N_Cols = 1) then
                    subA.[0, 0]
                else
                    let mutable d: float = 0.0
                    for j in 0 .. subA.N_Cols - 1 do
                        d <- d + subA.[0, j] * Math.Pow(-1.0, 2.0 + float j) * (AdvancedOps.SquareSubMatrix subA 0 j |> det)
                    d
        det A


    /// <summary>
    ///     This function computes the Euclidean norm of a Vector. This has been implemented
    ///     in Project A and is provided here for convenience
    /// </summary>
    /// <param name="v">
    ///    A Vector
    /// </param>
    /// <returns>
    ///     Euclidean norm, i.e. (\sum v[i]^2)^0.5
    /// </returns>
    static member VectorNorm (v : Vector) =
        let mutable n2 = 0.0
        for i in 0..v.Size-1 do
            n2 <- n2 + v.[i] * v.[i]
        sqrt n2
    

    /// <summary>
    ///     This function copies Vector 'v' as a column of matrix 'A'
    ///     at column position j.
    /// </summary>
    /// <param name="A">
    ///     An M-by-N matrix.
    /// </param>
    /// <param name="v">
    ///     Vector objects that must be copied in A.
    /// </param>
    /// <param name="j">
    ///     column number.
    /// </param>
    /// <returns>
    ///     An M-by-N matrix after modification.
    /// </returns>
    /// <exception cref="ArgumentException"></exception>
    static member SetColumn (A : Matrix) (v : Vector) (j : int) =
       //M-by-N matrix
       for i in 0 .. v.Size - 1 do
           A[i, j] <- v[i]
       A
       
    /// <summary>
    ///     This function computes the Gram-Schmidt process on a given matrix.
    /// </summary>
    /// <remarks>
    ///     See page 229 in "Linear Algebra for Engineers and Scientists"
    ///     by K. Hardy.
    /// </remarks>
    /// <param name="A">
    ///     An M-by-N matrix. All columns are implicitly assumed linear
    ///     independent.
    /// </param>
    /// <returns>
    ///     A tuple (Q,R) where Q is a M-by-N orthonormal matrix and R is an
    ///     N-by-N upper triangular matrix.
    /// </returns>
    static member GramSchmidt (A : Matrix) : Matrix * Matrix =
        let m = A.M_Rows
        let n = A.N_Cols
        let Q = Matrix(m, n)
        let R = Matrix(n, n)

        for j = 0 to n - 1 do
            let mutable v = A.Column j

            for i = 0 to j - 1 do
                R.[i, j] <- Q.Column i * v
                v <- v - Q.Column i * R.[i, j]
            
            R.[j, j] <- AdvancedOps.VectorNorm v

            let norm = R.[j, j]
            for k = 0 to m - 1 do
                Q.[k, j] <- v.[k] / norm

        (Q, R)
    
end