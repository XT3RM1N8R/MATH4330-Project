from basicFunctions import *

def BuildVandermonde(vector):
    """Builds a Vandermonde matrix for the given vector to the given 4th degree.

    Builds a vandermonde matrix for the given vector to the 4th degree, composed of quadratic polynomial power series for each element in the given vector.

    Args:
        vector: a vector that represents the independent set; the domain.

    Returns:
        A vandermonde matrix of the given vector to the 4th degree.
    """
    result = []
    for element in vector:
        temp = []
        for exponent in range(5): # Hardcoded for 4th degree
            temp.append(element**exponent)
        result.append(temp)
    return result

def ModifiedGramSchmidt(matrix):
    """Takes a given set of vectors and returns Q and R of QR Factorization

    Takes a given set of vectors and builds a temporary set, which is used in QR Factorization to build Q and R

    Args:
        matrix: A linearly independent set of vectors
    Returns:
        Q: a unitary matrix
        R: an upper triangular matrix
    """

    rowBounds = range(len(matrix))
    columnBounds = range(len(matrix[0]))

    V = [[0 for j in columnBounds] for i in rowBounds]
    for row in rowBounds:
        V[row] = matrix[row]
    V = Transpose(V)

    R = [[0 for j in columnBounds] for i in rowBounds]
    R = Transpose(R)
    Q = [[0 for j in columnBounds] for i in rowBounds]
    Q = Transpose(Q)
    for row in rowBounds:
        R[row][row] = TwoNorm(V[row])
        Q[row] = VectorScalarMultiplication(V[row], 1 / R[row][row])

        for k in range(row + 1, len(matrix)):
            R[k][row] = Dot(Q[row], V[k])
            temp = VectorScalarMultiplication(Q[row], R[k][row])
            V[k] = VectorSubtraction(V[k], temp)

    R = Transpose(R)
    Q = Transpose(Q)

    return [Q, R]

def BackSubstitution(matrix, vector):
    """Solves for the coefficients using back-substitution.

    Solves for the coefficients of the system by using back-substitution. Iterates the matrix in reverse.

    Args:
        matrix: an upper triangular matrix representing different partial equations for the system.
        vector: a vector representing the solution to the system.

    Returns:
        A vector of coefficients for the system.
    """
    result = vector
    for iterator in range(len(matrix[0])):
        x = (len(matrix[0]) - 1)
        summation = 0
        for k in range((x - iterator + 1), len(matrix)):
            summation += matrix[x - iterator][k] * result[k]
        result[x - iterator] = (vector[x - iterator] - summation) / matrix[x - iterator][x - iterator]
    return result


def Degree4Interpolation(coefficients, x):
    """Approximates a polynomial in the 4th degree.

    Get the value of the approximated polynomial for the given x.

    Args:
        coefficients: a vector of coefficients for the approximated polynomial.
        x: an arbitrary scalar that represents a value in the domain of the function.

    Returns:
        The value of the approximated polynomial.
    """
    result = (coefficients[0] +
              coefficients[1] * x +
              coefficients[2] * (x ** 2) +
              coefficients[3] * (x ** 3) +
              coefficients[4] * (x ** 4))
    return result


independentSet = [1, 2, 3, 4, 5]
dependentSet = [1j, 4j, 9j, 16j, 25j]

vandermondeMatrix = BuildVandermonde(independentSet)
[matrixQ, matrixR] = ModifiedGramSchmidt(vandermondeMatrix)

systemSolution = MatrixVectorMultiplication(ConjugateTranspose(matrixQ), dependentSet)
polynomialCoefficients = BackSubstitution(matrixR, systemSolution)

print(Degree4Interpolation(polynomialCoefficients, 5))

