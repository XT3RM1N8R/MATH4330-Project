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