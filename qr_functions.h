#ifndef QR_FUNCTIONS_H
#define QR_FUNCTIONS_H

#include <math.h>


int QR_decomposition(int n, double *matrix, double *result, double *d); // D - additional
void set_vector(int n, double *matrix, double *d, int k); // k - iteration
void multiplicate(int n, double *matrix, double *d, int k); // left multiplicate matrix on unitar
void build_result(int n, double *result, double *matrix, double *d);


void set_vector(int n, double *matrix, double *d, int k) {
    double s_k = 0;
    double norma;
    double norma_x;

    for (int i = k; i < n; i++)
        s_k += matrix[i * n + (k - 1)] * matrix[i * n + (k - 1)];


    norma = sqrt(
            matrix[(k - 1) * n + (k - 1)] * matrix[(k - 1) * n + (k - 1)] + s_k
    );

    matrix[(k - 1) * n + (k - 1)] -= norma;

    norma_x = sqrt(
            matrix[(k - 1) * n + (k - 1)] * matrix[(k - 1) * n + (k - 1)] + s_k
    );

    d[k - 1] = norma;

    for (int i = k - 1; i < n; i++)
        matrix[i * n + (k - 1)] /= norma_x;
}


void multiplicate(int n, double *matrix, double *d, int k) {
    for (int j = k; j < n; j++) {
        double s_p = 0;

        for (int i = k - 1; i < n; i++)
            s_p += matrix[i * n + j] * matrix[i * n + (k - 1)];
        s_p *= 2;

        for (int i = k - 1; i < n; i++)
            matrix[i * n + j] -= s_p * matrix[i * n + (k - 1)];
    }
}


void build_result(int n, double *result, double *matrix, double *d) {
    for (int k = 1; k < n; k++) {
        for (int j = 0; j < n; j++) {
            double s_p = 0;

            for (int i = k - 1; i < n; i++)
                s_p += result[i * n + j] * matrix[i * n + (k - 1)];
            s_p *= 2;

            for (int i = k - 1; i < n; i++)
                result[i * n + j] -= s_p * matrix[i * n + (k - 1)];
        }
    }


    for (int s = 0; s < n; s ++)
        matrix[s * n + s] = d[s];

    for (int j = 0; j < n - 1; j ++)
        for (int i = j + 1; i < n; i++)
            matrix[i * n + j] = 0.0;

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++)
            d[i] = result[i * n + j];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int l = 0; l < n; l++)
                sum += d[l] * matrix[i * n + l];
            result[i * n + j] = sum;
        }
    }
}


int QR_decomposition(int n, double *matrix, double *result, double *d) {
    int check_triangle = 1;
    for (int k = 1; k < n; k++) {
        int is_already = 1;
        for (int i = k; i < n; i++)
            if (matrix[i * n + (k - 1)] != 0) {
                is_already = 0;
                check_triangle = 0;
            }
        if (is_already == 1) {
            d[k - 1] = matrix[(k - 1) * n + (k - 1)];
            continue;
        }
        set_vector(n, matrix, d, k);
        multiplicate(n, matrix, d, k);
    }
    d[n - 1] = matrix[(n - 1) * n + (n - 1)];
    printf("R\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j > i)
                printf("%lf ", matrix[i * n + j]);
            if (j == i)
                printf("%lf ", d[i]);
            if (j < i)
                printf("0.0 ");
        }
        printf("\n");
    }
    printf("\n\n");

    for (int i = 0; i < n; i++)
        if (fabs(d[i]) < exp(-5))
            return 1;

    // finding R^-1
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++)
            matrix[i * n + j] /= d[i];
        result[i * n + i] /= d[i];
    }
    result[n * n - 1] /= d[n - 1];

    for (int j = n - 1; j > 0; j--)
        for (int i = j - 1; i >= 0; i--)
            for (int s = n - 1; s >= j; s--)
                result[i * n + s] -= result[j * n + s] * matrix[i * n + j];

    if (check_triangle == 1)
        return 0;

    printf("R^-1\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%lf ", result[i * n + j]);
        printf("\n");
    }
    printf("\n\n");
    // writing R^-1 to the matrix
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++) {
            matrix[i * n + j] = result[i * n + j];
            result[i * n + j] = 0.0;
        }

    for (int i = 0; i < n; i++) {
        d[i] = result[i * n + i];
        result[i * n + i] = 1.0;
    }

    // now R is a E-matrix, we can build Q matrix
    build_result(n, result, matrix, d);

    return 0;

}

#endif //QR_FUNCTIONS_H
