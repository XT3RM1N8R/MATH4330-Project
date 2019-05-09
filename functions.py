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
        for exponent in range(5):
            temp.append(element**exponent)
        result.append(temp)
    return result

def Transpose(matrix):
    """Get the transpose of a matrix.

    Get the transpose of a matrix or vector by rotating the column and row elements around the diagonal.

    Args:
        matrix: a matrix or vector.

    Returns:
        The transpose of the matrix.
    """
    result = []
    for iterator in range(len(matrix[0])):
        temp = []
        for element in range(len(matrix)):
            temp.append(matrix[element][iterator])
        result.append(temp)
    return result

def Conjugate(scalar):
    """Get the conjugate of a scalar.

    Get the conjugate of a scalar by subtracting the imaginary part from the real part.

    Args:
        scalar: a scalar.

    Returns:
        The conjugate of the scalar.
    """
    result = scalar.real
    result = result - scalar.imag
    return result

def ConjugateMatrix(matrix):
    """Get the conjugate of the matrix.

    Get the conjugate of the matrix by getting the conjugate of each element of the matrix.

    Args:
        matrix: a matrix.

    Returns:
        The conjugate of the matrix.
    """
    result = matrix
    for row in result:
        for element in row:
            Conjugate(element)
    return result

def ConjugateTranspose(matrix):
    """Gets the conjugate transpose of a matrix.

    Gets the conjugate transpose of a matrix by getting the transpose and then finding the conjugate of the matrix.

    Args:
        matrix: a matrix.

    Returns:
        The conjugate trasnpose of the matrix.
    """
    result = Transpose(matrix)
    for iterator in range(len(matrix[0])):
        for element in range(0, len(matrix)):
            result[iterator][element] = Conjugate(result[iterator][element])
    return result

def Dot(vector1, vector2):
    """Computes the scalar dot product of two vectors.

    This function computes a scalar which is the sum of all of the products between each pair of elements across both vectors.
    Args:
      vector1: A list of real numbers representing a vector.
      vector2: A list of real numbers of the same dimensions as vector1 also representing a vector.

    Returns:
      The dot product of the inputs.
    """

    result = 0
    for iterator in range(len(vector1)):
        result += (vector1[iterator] * vector2[iterator])
    return result

def TwoNorm(vector):
  """Calculates the 2-norm of a vector.

  Sums the squares of the elements of a given vector and returns the root of the sum.

  Args:
    vector: A list of numbers representing a vector.
  Returns:
    A scalar which is the 2-norm of the given vector.
  """
  result = 0
  for element in range(len(vector)):
    result = result + (vector[element]**2)
  result = result**(1/2)
  return result

def VectorScalarMultiplication(vector, scalar):
  """Computes the vector product between a scalar and each element of the vector.

  This function computes a vector that is the result of multiplying the given scalar by each element in the vector.

  Args:
    scalar: A single number
    vector: A list of numbers representing a vector.

  Returns:
    A list of numbers representing the desired scalar vector multiplication.
  """

  result = [0] * len(vector)
  for iterator in range(len(result)):
    result[iterator] = scalar*vector[iterator]
  return result

def MatrixVectorMultiplication(matrix, vector):
  """Computes the vector product between a matrix and a vector.

  This function computes a vector that is the result of computing the dot product of each row vector in the matrix by the given vector.

  Args:
    matrix: A list of vectors
    vector: A list of numbers representing a vector
  Returns:
    A vector that is the parallel dot product between each row vector in the matrix and the given vector.
  """
  
  result = [0] * len(matrix)
  for element in range(len(matrix)):
    result[element] = Dot(vector, matrix[element])
  return result

def vectorAdd(vector1, vector2):
  """Computes the vector sum of two vectors.

  This function computes the vector sum of two vectors by adding the respective elements of each vector together.

  Args:
    vector1: A list of real numbers representing a vector.
    vector2: A list of real numbers of the same dimensions as vector 1 also representing a vector.
  Returns:
    A list of real numbers representing the vector sum of the two given vectors.
  """

  result = [0]*len(vector1)
  for iterator in range(len(vector1)):
    result[iterator] = vector1[iterator] + vector2[iterator]
  return result

def VectorSubtraction(vector1, vector2):
  """Computes the vector sum of two vectors.

  This function computes the vector sum of two vectors by adding the respective elements of each vector together.

  Args:
    vector1: A list of real numbers representing a vector.
    vector2: A list of real numbers of the same dimensions as vector 1 also representing a vector.
  Returns:
    A list of real numbers representing the vector sum of the two given vectors.
  """

  result = [0] * len(vector1)
  for iterator in range(len(vector1)):
    result[iterator] = vector1[iterator] - vector2[iterator]
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

        for k in range(row+1, len(matrix)):
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
        x = (len(matrix[0])-1)
        summation = 0
        for k in range((x-iterator+1), len(matrix)):
            summation += matrix[x-iterator][k] * result[k]
        result[x-iterator] = (vector[x-iterator] - summation) / matrix[x-iterator][x-iterator]
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
              coefficients[2] * (x**2) +
              coefficients[3] * (x**3) +
              coefficients[4] * (x**4))
    return result

independentSet = [1j, 2j, 3j, 4j, 5j]
dependentSet = [1j, 4j, 9j, 16j, 25j]
degree = 4

vandermondeMatrix = BuildVandermonde(independentSet)
[matrixQ, matrixR] = ModifiedGramSchmidt(vandermondeMatrix)

systemSolution = MatrixVectorMultiplication(ConjugateTranspose(matrixQ), dependentSet)
polynomialCoefficients = BackSubstitution(matrixR, systemSolution)

print(Degree4Interpolation(polynomialCoefficients, 1))

