#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

#define EPSILON 1e-4
#define MAX_ITER 300
#define BETA 0.5


double** allocate_matrix(int n, int m) {
    int i;
    double** matrix = (double**)malloc(n * sizeof(double*));
    if (matrix == NULL) {
        return NULL;
    }
    for (i = 0; i < n; i++) {
        matrix[i] = (double*)calloc(m, sizeof(double));
        if (matrix[i] == NULL) {
            while (--i >= 0) {
                free(matrix[i]);
            }
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}


void free_matrix(double** matrix, int n) {
    int i;
    if (matrix == NULL) {
        return;
    }
    for (i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}


double** copy_matrix(double** src, int n, int m) {
    int i, j;
    double** dst = allocate_matrix(n, m);
    if (dst == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            dst[i][j] = src[i][j];
        }
    }
    return dst;
}


double euclidean_dist_squared(double* point1, double* point2, int dim) {
    double dist = 0.0;
    int i;
    for (i = 0; i < dim; i++) {
        dist += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    }
    return dist;
}


double** calculate_similarity_matrix(double** data, int n, int dim) {
    double** A = allocate_matrix(n, n);
    int i, j;
    double dist;
    if (A == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                dist = euclidean_dist_squared(data[i], data[j], dim);
                A[i][j] = exp(-dist / 2);
            }
        }
    }
    return A;
}


double** calculate_diagonal_degree_matrix(double** similarity, int n) {
    double** D = allocate_matrix(n, n);
    int i, j;
    double sum;
    if (D == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        sum = 0.0;
        for (j = 0; j < n; j++) {
            sum += similarity[i][j];
        }
        D[i][i] = sum;
    }
    return D;
}


double** calculate_normalized_similarity(double** A, int n) {
    double** D = calculate_diagonal_degree_matrix(A, n);
    double** W = allocate_matrix(n, n);
    int i, j;
    double di, dj;
    if (W == NULL) {
        free_matrix(D, n);
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            di = D[i][i];
            dj = D[j][j];
            if (di == 0.0 || dj == 0.0) {
                W[i][j] = 0.0;
            } else {
                W[i][j] = A[i][j] / (sqrt(di) * sqrt(dj));
            }
        }
    }
    free_matrix(D, n);
    return W;
}



double frobenius_norm_dist(double** A, double** B, int n, int m) {
    int i, j;
    double diff;
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            diff = A[i][j] - B[i][j];
            sum += diff * diff;
        }
    }
    return sum;
}


double** matrix_mul(double** A, double** B, int n, int m, int l) {
    int i, j, k;
    double** C = allocate_matrix(n, l);
    if (C == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < l; j++) {
            for (k = 0; k < m; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}


double** matrix_transpose(double** A, int n, int m) {
    int i, j;
    double** B = allocate_matrix(m, n);
    if (B == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            B[j][i] = A[i][j];
        }
    }
    return B;
}


double** symnmf(double** H, double** W, int n, int k) {
    int t, i, j;
    double numerator, denominator, change;
    double **WH, **Ht, **HHt, **HHtH;
    double** H_next = copy_matrix(H, n, k);

    for (t = 0; t < MAX_ITER; t++) {
        WH   = matrix_mul(W,  H,   n, n, k);
        Ht   = matrix_transpose(H, n, k);
        HHt  = matrix_mul(H,  Ht,  n, k, n);
        HHtH = matrix_mul(HHt, H,  n, n, k);

        for (i = 0; i < n; i++)
          for (j = 0; j < k; j++) {
            numerator   = WH[i][j];
            denominator = HHtH[i][j];
            if (denominator == 0) { printf("An Error Has Occurred\n"); exit(1);}
            H_next[i][j] = H[i][j] * (1 - BETA + BETA * (numerator / denominator));
          }

        free_matrix(WH, n);    free_matrix(Ht, k);
        free_matrix(HHt, n);   free_matrix(HHtH, n);
        
        change = frobenius_norm_dist(H_next, H, n, k);

        for (i = 0; i < n; i++)
          for (j = 0; j < k; j++)
            H[i][j] = H_next[i][j];

        if (change < EPSILON)
            break;
    }

    free_matrix(H_next, n);
    return H;
}



int main(int argc, char* argv[]) {
    char *goal, *file_name, *token, buffer[1024];
    FILE *fp;
    int i, j, n = 0, d = 0, lines = 0;
    double **data, **result;
    if (argc != 3) {
        printf("An Error Has Occurred\n");
        return 1; }
    goal = argv[1];
    file_name = argv[2];
    fp = fopen(file_name, "r");
    if (fp == NULL) {
        printf("An Error Has Occurred\n");
        return 1; }
    while (fgets(buffer, sizeof(buffer), fp)) {
        if (lines == 0) {
            token = strtok(buffer, ",");
            while (token != NULL) {
                d++;
                token = strtok(NULL, ","); }}
        lines++; }
    n = lines;
    rewind(fp);
    data = allocate_matrix(n, d);
    if (data == NULL) {
        printf("An Error Has Occurred\n");
        fclose(fp);
        return 1; }
    for (i = 0; i < n; i++) {
        fgets(buffer, sizeof(buffer), fp);
        token = strtok(buffer, ",");
        j = 0;
        while (token != NULL && j < d) {
            data[i][j] = atof(token);
            token = strtok(NULL, ",");
            j++; }}
    fclose(fp);
    if (strcmp(goal, "sym") == 0) {
        result = calculate_similarity_matrix(data, n, d);
    } else if (strcmp(goal, "ddg") == 0) {
        double **sim = calculate_similarity_matrix(data, n, d);
        result = calculate_diagonal_degree_matrix(sim, n);
        free_matrix(sim, n);
    } else if (strcmp(goal, "norm") == 0) {
        double **sim = calculate_similarity_matrix(data, n, d);
        result = calculate_normalized_similarity(sim, n);
        free_matrix(sim, n);
    } else {
        printf("An Error Has Occurred\n");
        free_matrix(data, n);
        return 1; }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%.4f", result[i][j]);
            if (j < n - 1) {
                printf(","); }}
        printf("\n"); }
    free_matrix(data, n);
    free_matrix(result, n);
    return 0; }