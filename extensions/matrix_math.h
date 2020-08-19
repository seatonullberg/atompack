#include <stdlib.h>

// Populates a matrix with the matrix multiplication product of two matrices.
void matmul(double a[][3], double b[][3], size_t length, double out[][3])
{
    float t;
    for (size_t i = 0; i < length; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            t = 0.0;
            for (int k = 0; k < 3; k++)
            {
                t += a[i][k] * b[k][j];
            }
            out[i][j] = t;
        }
    }
}

// Populates an array with the transpose of an array.
void transpose(double a[][3], size_t length, double out[][3])
{
    for (size_t i = 0; i < length; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            out[i][j] = a[j][i];
        }
    }
}