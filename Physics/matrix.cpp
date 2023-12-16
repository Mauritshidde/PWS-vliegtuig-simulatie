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
