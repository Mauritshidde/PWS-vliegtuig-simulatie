#include "matrix.h"

std::vector<std::vector<float>> invertMatrix(std::vector<std::vector<float>> matrix, int size)
{

    float a[size * size][size * size], x[size * size], ratio;
    int i, j, k, n;
    n = size;

    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {
            a[i][j] = matrix.at(i - 1).at(j - 1);
        }
    }

    /* Augmenting Identity Matrix of Order n */
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {
            if (i == j)
            {
                a[i][j + n] = 1;
            }
            else
            {
                a[i][j + n] = 0;
            }
        }
    }

    /* Applying Gauss Jordan Elimination */
    for (i = 1; i <= n; i++)
    {
        if (a[i][i] == 0.0)
        {
            std::cout << "Mathematical Error!";
            exit(0);
        }
        for (j = 1; j <= n; j++)
        {
            if (i != j)
            {
                ratio = a[j][i] / a[i][i];
                for (k = 1; k <= 2 * n; k++)
                {
                    a[j][k] = a[j][k] - ratio * a[i][k];
                }
            }
        }
    }
    /* Row Operation to Make Principal Diagonal to 1 */
    for (i = 1; i <= n; i++)
    {
        for (j = n + 1; j <= 2 * n; j++)
        {
            a[i][j] = a[i][j] / a[i][i];
        }
    }
    /* Displaying Inverse Matrix */
    std::vector<std::vector<float>> invertedMatrix = matrix;
    std::cout << std::endl
              << "Inverse Matrix is:" << std::endl;

    int i2 = 0;
    int j2 = 0;

    for (i = 1; i <= n; i++)
    {
        for (j = n + 1; j <= 2 * n; j++)
        {
            invertedMatrix.at(i - 1).at(j - n - 1) = a[i][j];
            std::cout << a[i][j] << "\t";
            j2++;
        }
        i2++;
        std::cout << std::endl;
    }

    return invertedMatrix;
}

// int main()
// {
//     std::vector<std::vector<float>> matrix = {
//         {1, 2, 3, 4, 5, 6}, {5, 3, 41, 2, 3, 4}, {3, 4, 2, 1, 2, 0}, {23, 4, 45, 3, 2, 1}, {12, 32, 3224, 4, 2, 42}, {42, 42, 42, 42, 2, 2}};
//     matrix = invertMatrix(matrix, matrix.size());

//     for (int i = 0; i < matrix.size(); i++)
//     {
//         for (int j = 0; j < matrix.size(); j++)
//         {
//             std::cout << matrix.at(i).at(j) << " ";
//         }
//         std::cout << std::endl;
//     }
//     return (0);
// }

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