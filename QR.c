#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <time.h>



void QR_decomposition(int n, double *matrix, double *result, double *d); // D - additional
void set_vector(int n, double *matrix, double *d, int k); // k - iteration
void multiplicate(int n, double *matrix, double *d, int k); // left multiplicate matrix on unitar
void build_result(int n, double *result, double *matrix, double *d);


int main(int argc, char **argv) {
	if (argc != 2) {
		printf("You should write filename\n");
		return -1;
	}

	FILE *fi = fopen(argv[1], "r");
	if (fi == NULL) {
		printf("incorrect filename\n");
		return -1;
	}

	int n;
	if (fscanf(fi, "%d", &n) != 1) {
		printf("Error in matrix size\n");
		return -1;
	}

	double *matrix = (double *)malloc(n * n * sizeof(double));
	double *result = (double *)malloc(n * n * sizeof(double));
	double *d = (double *)malloc(n * sizeof(double));

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

	QR_decomposition(n, matrix, result, d);

	time_t end = clock();

	printf("execution time: %f sec\n", (double)(end - start) / CLOCKS_PER_SEC);
	free(matrix);
	free(result);
	free(d);
	return 0;
}


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
	for (int k = 1; k <= n; k++) {
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

	printf("multi matrix:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			printf("%lf ", result[i * n + j]);
		printf("\n");
	}
	printf("\n");


}


void QR_decomposition(int n, double *matrix, double *result, double *d) {
	for (int k = 1; k <= n; k++) {
		set_vector(n, matrix, d, k);
		multiplicate(n, matrix, d, k);
	}

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

}


