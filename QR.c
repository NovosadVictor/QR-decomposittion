#include <malloc.h>
#include <time.h>

#include "qr_functions.h"
#include "formulas.h"



int main(int argc, char **argv) {
	if (argc != 3) {
		printf("You should choose mode(1 - formula, 2 - file)\n");
        printf("And write matrix size(for mode 1) of filename(for mode 2)\n");
		return -1;
	}

    if (atoi(argv[1]) == 2) {
        FILE *fi = fopen(argv[2], "r");
        if (fi == NULL) {
            printf("incorrect filename\n");
            return -1;
        }

        int n;
        if (fscanf(fi, "%d", &n) != 1) {
            printf("Error in matrix size\n");
            return -1;
        }

        double *matrix = (double *) malloc(n * n * sizeof(double));
        double *result = (double *) malloc(n * n * sizeof(double));
        double *d = (double *) malloc(n * sizeof(double));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double elem;
                if (fscanf(fi, "%lf", &elem) != 1) {
                    printf("Incorrect matrix in file\n");
                    return -1;
                }
                matrix[i * n + j] = elem;
            }
        char c;
        if (fscanf(fi, "%c", &c) > 0) {
            printf("NOT EOF\n");
            return -1;
        }

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    result[i * n + j] = 1.0;
                else
                    result[i * n + j] = 0.0;
            }

        time_t start = clock();

        if (QR_decomposition(n, matrix, result, d) == 1) {
            printf("\nDegenerate matrix\n");
            return -1;
        }

        time_t end = clock();
        print(result, n);

        printf("execution time: %f sec\n", (double) (end - start) / CLOCKS_PER_SEC);
        free(matrix);
        free(result);
        free(d);
        fclose(fi);
    }
    if (atoi(argv[1]) == 1) {
        int n = atoi(argv[2]);

        double *matrix = (double *) malloc(n * n * sizeof(double));
        double *result = (double *) malloc(n * n * sizeof(double));
        double *d = (double *) malloc(n * sizeof(double));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                matrix[i * n + j] = formula(i, j);


        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    result[i * n + j] = 1.0;
                else
                    result[i * n + j] = 0.0;
            }

        time_t start = clock();

        if (QR_decomposition(n, matrix, result, d) == 1) {
                printf("\nDegenerate matrix\n");
                return -1;
            }

        time_t end = clock();
        print(result, n);

        printf("execution time: %f sec\n", (double) (end - start) / CLOCKS_PER_SEC);
        free(matrix);
        free(result);
        free(d);
    }
	return 0;
}

