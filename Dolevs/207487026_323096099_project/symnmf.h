#ifndef SYMNMF_H
#define SYMNMF_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>


/* -----------------------------------------------------------------------------
 * Macros
 * -------------------------------------------------------------------------- */

#define ERROR_MESSAGE "An Error Has Occurred"
#define E 2.718281828459045

/* -----------------------------------------------------------------------------
 * Struct Definitions
 * -------------------------------------------------------------------------- */

typedef struct cord {
    double value;           /* Coordinate value */
    struct cord *next;      /* Pointer to the next coordinate */
} cord;

typedef struct vector {
    struct vector *next;    /* Pointer to the next vector */
    cord *cords;            /* Pointer to the list of coordinates */
} vector;

typedef struct matrix {
    double **mat;           /* Pointer to the matrix data */
    int rows;               /* Number of rows */
    int cols;               /* Number of columns */
} matrix;

/* -----------------------------------------------------------------------------
 * Memory Management and Utilities
 * -------------------------------------------------------------------------- */

/** Terminates the program with a given error message. */
void terminate(const char *message);

/** Allocates and returns a new vector node. Returns NULL on failure. */
vector *create_vector_node(void);

/** Allocates and returns a new coordinate node. Returns NULL on failure. */
cord *create_cord_node(void);

/** Frees a linked list of coordinates. */
void free_cords(cord *head_cord);

/** Frees a linked list of vectors (and their internal cords). */
void free_vectors(vector *head_vec);

/** Frees both vectors and cords. */
void cleanup_all(vector *vec, cord *cord);

/** Frees a dynamically allocated matrix (as array of arrays). */
void free_matrix(double **matrix, int rows);

/** Frees a matrix inside a matrix struct. */
void free_matrix_ptr(matrix *curr);

/** Frees multiple matrices passed as variable arguments. */
void free_matrixes(matrix *first, ...);

/* -----------------------------------------------------------------------------
 * Matrix Initialization
 * -------------------------------------------------------------------------- */

/** Initializes and returns a zero matrix (array of arrays). Returns NULL on failure. */
double **init_zeros_mat(int rows, int cols);

/** Fills a matrix from linked list structure. Terminates on failure. */
double **fill_matrix(vector *head_vec, int rows, int cols);

/** Returns an empty matrix struct (rows, cols = 0, mat = NULL). */
matrix empty_matrix(void);

/** Clones a matrix (deep copy). Returns empty matrix on failure. */
matrix clone_matrix(matrix *curr);

/** Initializes a matrix struct with a zero matrix of given size. */
matrix init_zero_struct_matrix(int rows, int cols);

/** Checks if matrix is valid (non-NULL and positive dimensions). */
int is_valid_matrix(matrix *curr);

/* -----------------------------------------------------------------------------
 * Matrix Operations
 * -------------------------------------------------------------------------- */

/** Multiplies two matrices. Returns empty matrix on failure. */
matrix matrix_mult(matrix mat1, matrix mat2);

/** Transposes a matrix. Returns empty matrix on failure. */
matrix transpose(matrix data);

/** Raises diagonal elements to the power of -0.5. Returns empty matrix on failure. */
matrix power_minus_half_diag(matrix data);

/** Subtracts right_mat from left_mat. Returns empty matrix on failure. */
matrix substruct_matrix(matrix left_mat, matrix right_mat);

/** Computes the squared Frobenius norm between two matrices. Returns -1 on failure. */
double squared_frobenius_norm(matrix left_mat, matrix right_mat);

/** Prints matrix to stdout. */
void print_matrix(matrix data);

/* -----------------------------------------------------------------------------
 * File Reading
 * -------------------------------------------------------------------------- */

/** Parses input file into vectors and cords. Returns 1 on success, 0 on failure. */
int parse_file_loop(FILE *fp, vector **curr_vec, cord **curr_cord, cord **head_cord, int *rows, int *cols);

/** Reads data from file and prints resulting matrix. Returns 0 on success. */
matrix read_data(const char *path);

/* -----------------------------------------------------------------------------
 * symNMF Core Functions
 * -------------------------------------------------------------------------- */

/** Computes the squared Euclidean distance between two rows. */
double squared_euqlidean_distance(matrix matrix_data, int row_index1, int row_index2);

/** Computes similarity matrix. Returns empty matrix on failure. */
matrix sym(matrix data);

/** Computes degree diagonal matrix. Returns empty matrix on failure. */
matrix ddg(matrix data);

/** Computes normalized matrix: D^(-0.5) * A * D^(-0.5). */
matrix norm(matrix data);

/** Performs one update iteration of H. Returns updated matrix or empty on failure. */
matrix update_H(matrix lap_mat, matrix curr, double beta);

/** Iteratively converges to updated H matrix. Returns final H or empty on failure. */
matrix convergence(matrix h_mat, matrix w_mat, double beta, double epsilon, int max_iter);

#endif /* SYMNMF_H */
