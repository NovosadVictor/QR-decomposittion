#ifndef QR_FORMULAS_H
#define QR_FORMULAS_H

#include <stdio.h>
#include <stdlib.h>


double formula(int i, int j) {
    return (i + j + 1) / 10 ? i + j > 3 : (i + 1) * (j - 7) / 5.0;
}

void print(double *matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
            printf("%lf ", matrix[i * size + j]);
        printf("\n");
    }
    printf("\n");
}

#endif //QR_FORMULAS_H
