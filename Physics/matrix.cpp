#include "matrix.h"

float getDeterminant(const std::vector<std::vector<float>> vect)
{
      if (vect.size() != vect[0].size())
      {
            throw std::runtime_error("Matrix is not quadratic");
      }
      int dimension = vect.size();

      if (dimension == 0)
      {
            return 1;
      }

      if (dimension == 1)
      {
            return vect[0][0];
      }

      // Formula for 2x2-matrix
      if (dimension == 2)
      {
            return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
      }

      float result = 0;
      int sign = 1;
      for (int i = 0; i < dimension; i++)
      {

            // Submatrix
            std::vector<std::vector<float>> subVect(dimension - 1, std::vector<float>(dimension - 1));
            for (int m = 1; m < dimension; m++)
            {
                  int z = 0;
                  for (int n = 0; n < dimension; n++)
                  {
                        if (n != i)
                        {
                              subVect[m - 1][z] = vect[m][n];
                              z++;
                        }
                  }
            }

            // recursive call
            result = result + sign * vect[0][i] * getDeterminant(subVect);
            sign = -sign;
      }

      return result;
}

std::vector<std::vector<float>> getTranspose(const std::vector<std::vector<float>> matrix1)
{

      // Transpose-matrix: height = width(matrix), width = height(matrix)
      std::vector<std::vector<float>> solution(matrix1[0].size(), std::vector<float>(matrix1.size()));

      // Filling solution-matrix
      for (size_t i = 0; i < matrix1.size(); i++)
      {
            for (size_t j = 0; j < matrix1[0].size(); j++)
            {
                  solution[j][i] = matrix1[i][j];
            }
      }
      return solution;
}

std::vector<std::vector<float>> getCofactor(const std::vector<std::vector<float>> vect)
{
      if (vect.size() != vect[0].size())
      {
            throw std::runtime_error("Matrix is not quadratic");
      }

      std::vector<std::vector<float>> solution(vect.size(), std::vector<float>(vect.size()));
      std::vector<std::vector<float>> subVect(vect.size() - 1, std::vector<float>(vect.size() - 1));

      for (std::size_t i = 0; i < vect.size(); i++)
      {
            for (std::size_t j = 0; j < vect[0].size(); j++)
            {

                  int p = 0;
                  for (size_t x = 0; x < vect.size(); x++)
                  {
                        if (x == i)
                        {
                              continue;
                        }
                        int q = 0;

                        for (size_t y = 0; y < vect.size(); y++)
                        {
                              if (y == j)
                              {
                                    continue;
                              }

                              subVect[p][q] = vect[x][y];
                              q++;
                        }
                        p++;
                  }
                  solution[i][j] = pow(-1, i + j) * getDeterminant(subVect);
            }
      }
      return solution;
}

std::vector<std::vector<float>> getInverse(const std::vector<std::vector<float>> vect)
{
      if (getDeterminant(vect) == 0)
      {
            throw std::runtime_error("Determinant is 0");
      }

      float d = 1.0 / getDeterminant(vect);
      std::vector<std::vector<float>> solution(vect.size(), std::vector<float>(vect.size()));

      for (size_t i = 0; i < vect.size(); i++)
      {
            for (size_t j = 0; j < vect.size(); j++)
            {
                  solution[i][j] = vect[i][j];
            }
      }

      solution = getTranspose(getCofactor(solution));

      for (size_t i = 0; i < vect.size(); i++)
      {
            for (size_t j = 0; j < vect.size(); j++)
            {
                  solution[i][j] *= d;
            }
      }

      return solution;
}

std::vector<std::vector<float>> zeros(int width, int height)
{
      std::vector<std::vector<float>> vector;

      for (int i = 0; i < width; i++)
      {
            std::vector<float> helper;
            for (int j = 0; j < height; j++)
            {
                  helper.push_back(0);
            }
            vector.push_back(helper);
      }

      return vector;
}

float det(std::vector<std::vector<float>> A)
{
}

float determinantOfMatrix(std::vector<std::vector<float>> mat, int n)
{
      float num1, num2, det = 1, index, total = 1;

      float temp[n + 1];

      for (int i = 0; i < n; i++)
      {
            index = i;
            while (index < n && mat.at(index).at(i) == 0)
            {
                  index++;
            }
            if (index == n)
            {
                  continue;
            }
            if (index != i)
            {
                  for (int j = 0; j < n; j++)
                  {
                        std::swap(mat.at(index).at(j), mat.at(i).at(j));
                  }
                  det = det * pow(-1, index - i);
            }

            for (int j = 0; j < n; j++)
            {
                  temp[j] = mat.at(i).at(j);
            }
            for (int j = i + 1; j < n; j++)
            {
                  num1 = temp[i];         // value of diagonal element
                  num2 = mat.at(j).at(i); // value of next row element

                  for (int k = 0; k < n; k++)
                  {
                        mat.at(j).at(k) = (num1 * mat.at(j).at(k)) - (num2 * temp[k]);
                  }
                  total = total * num1; // Det(kA)=kDet(A);
            }
      }

      for (int i = 0; i < n; i++)
      {
            det = det * mat.at(i).at(i);
      }
      return (det / total); // Det(kA)/k=Det(A);
}

std::vector<float> linspace(int startX, int endX, int steps)
{
      float stepSize = (endX - startX) / (steps - 1);
      std::vector<float> coords;

      for (int i = 0; i < steps; i++)
      {
            coords.push_back(startX + (stepSize * i));
      }

      return coords;
}