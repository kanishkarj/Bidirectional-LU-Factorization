# Bidirectional LU Factorization


### Prepared For
    Dr. Surya Prakash


### Prepared By
    Pushpendra Kumar : 160001046
    Kanishkar J : 160001028


## Introduction


    In numerical analysis and linear algebra, lower–upper (LU) decomposition or factorization factors a matrix as the product of a lower triangular matrix and an upper triangular matrix. The product sometimes includes a permutation matrix as well. The LU decomposition can be viewed as the matrix form of Gaussian elimination. Computers usually solve square systems of linear equations using the LU decomposition, and it is also a key step when inverting a matrix or computing the determinant of a matrix.


    As, we can see how important LU decomposition is in scientific computations. And the computing power of a single processor is not always enough for processing very big data. Hence, In this project we shall try to Parallelize the LU decomposition method (Bidirectional LU factorization), and we shall also compare its performance metrics with the normal LU decomposition method. 


    In this project we will implement a parallel _bidirectional algorithm_, based on _LU factorization_, for the solution of general sparse system of linear equations having _non symmetric_ coefficient matrix. As with the sparse symmetric systems, the numerical factorization phase of our algorithm is carried out in such a manner that the entire back substitution component of the substitution phase is replaced by a single step division. However, due to absence of symmetry, important differences arise in the ordering technique, the symbolic factorization phase, and message passing during numerical factorization phase. The bidirectional substitution phase for solving general sparse systems is the same as that for sparse symmetric systems. 


## Project Roadmap


### Implementation 

First, We will be implementing the sequential and parallel algorithms for LU Decomposition problem in c++ using OpenMP/CUDA .  


### Verification 


    We shall then verify the outputs of both the procedures for different matrices and compare them.


### Analysis 



*   Compare Speed-ups for Matrices of different sizes and compare plot a graph of the data obtained.
*   Compare Speed-ups while increasing the number of processors used and plot a graph with it.


### Testing

We will test our algorithm on randomly generated matrices of various densities.


### Solve System of Linear Equations

We will use LU decomposition to solve a system of linear equations. Here we will be using both the sequential and the parallel algorithm and then note down the speed-up.


## Methodology


### LU Factorization 

It is always possible to factor a square matrix into a lower triangular matrix and an upper triangular matrix. That is, [A] = [L][U]. Doolittle's method provides an alternative way to factor A into an LU decomposition without going through the hassle of Gaussian Elimination. For a general n×n matrix A, we assume that an LU decomposition exists, and write the form of L and U explicitly. We then systematically solve for the entries in L and U from the equations that result from the multiplications necessary for A=LU.


``` cpp
int lower[n][n]
int upper[n][n];

// Decomposing matrix into Upper and Lower triangular matrix
for (int i = 0; i < n; i++) {

   // Upper Triangular
   for (int k = i; k < n; k++) {
       // Summation of L(i, j) * U(j, k)
       int sum = 0;
       for (int j = 0; j < i; j++)
           sum += (lower[i][j] * upper[j][k]);
       // Evaluating U(i, k)
       upper[i][k] = mat[i][k] - sum;
   }

   // Lower Triangular
   for (int k = i; k < n; k++) {
       if (i == k)
           lower[i][i] = 1; // Diagonal as 1
       else {
           // Summation of L(k, j) * U(j, i)
           int sum = 0;
           for (int j = 0; j < i; j++)
               sum += (lower[k][j] * upper[j][i]);
           // Evaluating L(k, i)
           lower[k][i] = (mat[k][i] - sum) / upper[i][i];
       }
   }
}
```


After implementing the sequential algorithm, We will parallelize the above mentioned algorithms to achieve better performance.


### Bidirectional Substitution



1.  Forward elimination: reduction to row echelon form. Using it one can tell whether there are no solutions, or unique solution, or infinitely many solutions.
2.  Back substitution: further reduction to reduced row echelon form.

We will be implementing methods for forward and backward substitution. Then we will be applying forward substitution and then backward substitution to solve the system of linear equations. 

