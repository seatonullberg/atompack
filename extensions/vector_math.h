#include <math.h>

// Populates an array with the sum of two arrays.
void add(double a[3], double b[3], double out[3])
{
    for (int i = 0; i < 3; i++)
    {
        out[i] = a[i] + b[i];
    }
}

// Populates an array with the cross product of two arrays.
void cross(double a[3], double b[3], double out[3])
{
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

// Populates an array with the quotient of an array and a scalar.
void div_scalar(double a[3], double b, double out[3])
{
    for (int i = 0; i < 3; i++)
    {
        out[i] = a[i] / b;
    }
}

// Returns the dot product of two arrays.
double dot(double a[3], double b[3])
{
    double prod = 0.0;
    for (int i = 0; i < 3; i++)
    {
        prod += a[i] * b[i];
    }
    return prod;
}

// Returns the euclidean distance between two arrays.
double euclidean(double a[3], double b[3])
{
    double sum = 0.0;
    for (int i = 0; i < 3; i++)
    {
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(sum);
}

// Returns the L2 norm of an array.
double norm(double a[3])
{
    double sum = 0.0;
    for (int i = 0; i < 3; i++)
    {
        sum += a[i] * a[i];
    }
    return sqrt(sum);
}

// Populates an array with the difference between two arrays.
void sub(double a[3], double b[3], double out[3])
{
    for (int i = 0; i < 3; i++)
    {
        out[i] = a[i] - b[i];
    }
}