#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include  <math.h>
#include <stdarg.h>
#include <string.h>
#include "symnmf.h"


#define ERROR_MESSAGE "An Error Has Occurred"
#define E 2.718281828459045


void terminate(const char *message)
{
    printf("%s\n", message);
    exit(1);
}

int terminate_main(const char *message)
{
    printf("%s\n", message);
    return -1;
}

/* -----------------------------------------------------------------------------
 * Vectors and Cords for Reading
 * ----------------------------------------------------------------------------- */
 
struct vector *create_vector_node(void)
{
    struct vector *v = NULL;
    v = malloc(sizeof(struct vector));
    if (v == NULL) {
        return NULL; }
    v->cords = NULL;
    v->next = NULL;
    return v;
}


struct cord *create_cord_node() {
    struct cord *cord = NULL;
    cord = malloc(sizeof(struct cord));
    if (cord == NULL)  {
        return NULL; }
    cord->next = NULL;
    cord->value = 0.0;
    return cord;
}

void free_cords(struct cord *head_cord) {
    struct cord *temp;

    while (head_cord != NULL) {
        temp = head_cord;
        head_cord = head_cord->next;
        free(temp);
    }
}

void free_vectors(struct vector *head_vec) {
    struct vector *temp;

    while (head_vec != NULL) {
        
        if (head_vec->cords != NULL) {
           free_cords(head_vec->cords);
        }
    
        temp = head_vec;
        head_vec = head_vec->next;
        free(temp);
    }
}

void cleanup_all(struct vector *vec, struct cord *cord) {
    free_vectors(vec);
    free_cords(cord);
}

/* -----------------------------------------------------------------------------
 * Matrix Functions
 * ----------------------------------------------------------------------------- */

void free_matrix(double **matrix, int rows) {
    int i;
    for (i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

double** init_zeros_mat(int rows, int cols) {
    int i, j;
    double **init_matrix = malloc(rows * sizeof(double *));
    if (init_matrix == NULL) {
        return NULL;  
    }
    for (i = 0; i < rows; i++) {
        init_matrix[i] = calloc(cols, sizeof(double));
        if (init_matrix[i] == NULL) {
            for (j = 0; j < i; j++) {
                free(init_matrix[j]);
            }
            free(init_matrix);
            return NULL;
        }
    }
    return init_matrix;
}

double** fill_matrix(struct vector *head_vec, int rows, int cols) {
    int i,j;
    struct vector *curr_vec = head_vec;
    struct cord *curr_cord;
    double **vectors_matrix;
    if(curr_vec == NULL || head_vec == NULL) {
        terminate(ERROR_MESSAGE);
    }
    vectors_matrix = malloc(rows * sizeof(double *));
    if (vectors_matrix == NULL) {
        terminate(ERROR_MESSAGE);
    }
    for (i = 0; i < rows; i++) {
        vectors_matrix[i] = malloc(cols * sizeof(double));
        if(vectors_matrix[i] == NULL) {
            for (j = 0; j < i; j++) {
                free(vectors_matrix[j]);
            }
            free(vectors_matrix);
            free(head_vec);
            terminate(ERROR_MESSAGE);
        }
        for (j = 0; j < cols; j++) {
            vectors_matrix[i][j] = 0.0;
        }
    }
    i = 0;
    while (curr_vec != NULL && curr_vec->cords != NULL && i < rows) {
        curr_cord = curr_vec->cords;
        j = 0;
        while (curr_cord != NULL && j < cols) {
            vectors_matrix[i][j] = curr_cord->value;
            curr_cord = curr_cord->next;
            j++;
        }
        i++;
        curr_vec = curr_vec->next;
    }
    return vectors_matrix;
}



int is_valid_matrix(matrix *curr) {
    return (curr != NULL && curr->rows > 0 && curr->cols > 0 && curr->mat != NULL);
}

void free_matrix_ptr(matrix *curr) {
    int i;

    if (!is_valid_matrix(curr)) {
        return;
    }

    for (i = 0; i < curr->rows; i++) {
        free(curr->mat[i]); 
    }

    free(curr->mat);
    curr->mat = NULL;
    curr->rows = 0;
    curr->cols = 0;
}

void free_matrix_ptr_old(matrix *curr) {
    int i;
    if ((curr == NULL) || (curr->mat == NULL))
    {
    return;
    }
    for (i = 0; i < curr->rows; i++) 
    {
        if (curr->mat[i] != NULL) 
        {
            free(curr->mat[i]);
        }
    }
    if (curr->mat != NULL)
    {
        free(curr->mat);
    }
    curr->mat = NULL;
}

void free_matrixes(matrix *first, ...) {
    va_list args;
    matrix *curr = first;

    va_start(args, first);

    while (curr != NULL) {
        free_matrix_ptr(curr);          
        curr = va_arg(args, matrix*);   
    }

    va_end(args);
}

matrix empty_matrix() {
    matrix m;
    m.mat = NULL;
    m.rows = 0;
    m.cols = 0;
    return m;
}

matrix clone_matrix(matrix *curr){   
    int i, j;
    double **data;
    matrix cloned;
    if ((curr == NULL) || (curr->mat == NULL)){
    return empty_matrix();
    }
    if (!is_valid_matrix(curr)){
        empty_matrix();
    }
    data = init_zeros_mat(curr->rows, curr->cols);
    if(data == NULL){
        return empty_matrix();
    }
    for(i=0; i<(curr->rows); i++)
    {
        for(j=0 ; j<(curr->cols); j++)
        {
            data[i][j] = curr->mat[i][j];
        }
    }
    cloned.mat = data;
    cloned.rows = curr->rows;
    cloned.cols = curr->cols;
    return cloned;
}
matrix init_zero_struct_matrix(int rows, int cols){
    matrix zero_matrix;
    zero_matrix.mat = init_zeros_mat(rows,cols);
    if (zero_matrix.mat == NULL){
        return empty_matrix();
    }
    zero_matrix.rows = rows;
    zero_matrix.cols = cols;
    return zero_matrix;
}

void print_matrix(matrix data) {
    int i, j;
    if (!is_valid_matrix(&data)) {
        return;
    }

    for (i = 0; i < data.rows; i++) {
        for (j = 0; j < data.cols; j++) {
            if (j == data.cols - 1) {
                printf("%.4f", data.mat[i][j]); 
            } else {
                printf("%.4f,", data.mat[i][j]);
            }
        }
        printf("\n");
    }
}


/* -----------------------------------------------------------------------------
 * Arithmetics
 * ----------------------------------------------------------------------------- */
matrix power_minus_half_diag(matrix data)
{   int i;
    matrix result;
    if (!is_valid_matrix(&data)) {
        return empty_matrix();
    }

    result = clone_matrix(&data);
    if (!is_valid_matrix(&result)){
        return empty_matrix();
    }
    for(i=0; i< result.rows; i++) {
        result.mat[i][i] = pow(result.mat[i][i], -0.5);
    }
    return result;

}

matrix matrix_mult(matrix mat1, matrix mat2){ 
    int i,j,k;
    double temp;
    matrix mult_matrix;
    if ((!is_valid_matrix(&mat1))||(!is_valid_matrix(&mat2)))
    {
        return empty_matrix();
    }
    if (mat1.cols != mat2.rows) {
        return empty_matrix();
    }
    mult_matrix = init_zero_struct_matrix(mat1.rows,mat2.cols);
    for (i=0; i<mult_matrix.rows;i++){
        for (j=0; j<mult_matrix.cols;j++){
            temp = 0;
            for (k=0 ; k<mat1.cols ; k++){
                temp += ((mat1.mat[i][k]) * (mat2.mat[k][j]));
            }
            mult_matrix.mat[i][j] = temp;
        }
    }
    return mult_matrix;
}

matrix transpose(matrix data) {
    int i, j;
    matrix result;

    if (!is_valid_matrix(&data)) {
        return empty_matrix();
    }

    result = init_zero_struct_matrix(data.cols, data.rows);
    if (!is_valid_matrix(&result)) {
        return empty_matrix();
    }

    for (i = 0; i < data.rows; i++) {
        for (j = 0; j < data.cols; j++) {
            result.mat[j][i] = data.mat[i][j];
        }
    }
    return result;
}



matrix substruct_matrix(matrix left_mat, matrix right_mat){
    int i, j;
    matrix result; 
    if ((left_mat.rows != right_mat.rows) || (left_mat.cols!= right_mat.cols) || !is_valid_matrix(&left_mat) || !is_valid_matrix(&right_mat)){
        return empty_matrix();
    }
    result = clone_matrix(&left_mat);
    if(!is_valid_matrix(&result)){
        return empty_matrix();
    }
    for (i = 0; i < left_mat.rows; i++){
        for (j = 0; j < left_mat.cols; j++){
            result.mat[i][j] = result.mat[i][j] - right_mat.mat[i][j];
        }
    }
    return result;

}

double squared_frobenius_norm(matrix left_mat, matrix right_mat){
    int i, j;
    double cnt = 0.0, diff;
    if ((left_mat.rows != right_mat.rows) || (left_mat.cols != right_mat.cols) || !is_valid_matrix(&left_mat) || !is_valid_matrix(&right_mat)) {
        return -1;
    }
    for (i = 0; i < left_mat.rows; i++) {
        for (j = 0; j < left_mat.cols; j++) {
            diff = left_mat.mat[i][j] - right_mat.mat[i][j];
            cnt += (diff * diff);
        }
    }
    return cnt; 
}

/* -----------------------------------------------------------------------------
 * Read Data
 * ----------------------------------------------------------------------------- */

 int parse_file_loop(FILE *fp, struct vector **curr_vec, struct cord **curr_cord, struct cord **head_cord, int *rows, int *cols) {
    double n;
    char c;
    struct vector *new_vec;
    int temp_cols = 0, first = 1;
    while (fscanf(fp, "%lf%c", &n, &c) == 2) {
        if (c == '\n') {
            (*curr_cord)->value = n;
            (*curr_vec)->cords = *head_cord;
            new_vec = create_vector_node();
            if (!new_vec) return 0;
            (*curr_vec)->next = new_vec;
            *curr_vec = new_vec;

            *head_cord = create_cord_node();
            if (!(*head_cord)) return 0;
            *curr_cord = *head_cord;

            (*rows)++;
            if (first) {
                *cols = temp_cols + 1;
                first = 0; 
            }
            temp_cols = 0; 
        } 
        else {
            temp_cols++;
            (*curr_cord)->value = n;
            (*curr_cord)->next = create_cord_node();
            if (!(*curr_cord)->next) return 0;
            *curr_cord = (*curr_cord)->next; 
        }
    }
    return 1;
}

matrix read_data(const char *path) {
    int rows = 0, cols = 0;
    double **result_matrix;
    matrix matrix_data;
    vector *head_vec, *curr_vec;
    cord *head_cord, *curr_cord;
    FILE *ifp = fopen(path, "r");
    if (!ifp){
        terminate(ERROR_MESSAGE);
    }
    head_vec = create_vector_node();
    head_cord = create_cord_node();
    curr_vec = head_vec;
    curr_cord = head_cord;
    if (!parse_file_loop(ifp, &curr_vec, &curr_cord, &head_cord, &rows, &cols)) {
        fclose(ifp);
        cleanup_all(head_vec, head_cord);
        terminate(ERROR_MESSAGE);
    }
    result_matrix = fill_matrix(head_vec, rows, cols);
    if (!result_matrix) {
        fclose(ifp);
        cleanup_all(head_vec, head_cord);
        terminate(ERROR_MESSAGE);
    }
    fclose(ifp);
    free_vectors(head_vec);
    free_cords(head_cord); 
    matrix_data.mat = result_matrix;
    matrix_data.rows = rows;
    matrix_data.cols = cols;
    if(!is_valid_matrix(&matrix_data)){
        terminate(ERROR_MESSAGE);
    }
    return matrix_data;
}

/* -----------------------------------------------------------------------------
 * 1.1
 * ----------------------------------------------------------------------------- */

double squared_euqlidean_distance(matrix matrix_data,int row_index1, int row_index2)
{
    double distance = 0, temp_num;
    int i;
    for(i =0 ; i < matrix_data.cols ; i ++){
        temp_num = ((matrix_data.mat[row_index1][i]) - (matrix_data.mat[row_index2][i]));
        distance += temp_num * temp_num;
    }
    return distance;
}

matrix sym(matrix data)
{
    matrix sym; 
    int i, j;
    double power = 0;
    if (!is_valid_matrix(&data)) {
        return empty_matrix();
    }
    sym = init_zero_struct_matrix(data.rows, data.rows);
    if (!is_valid_matrix(&sym)) {
        return empty_matrix();
    }

    for (i = 0; i < sym.rows; i++) {
        for (j = i + 1; j < sym.cols; j++) {
            power = -((squared_euqlidean_distance(data, i, j)) / 2);
            sym.mat[i][j] = exp(power); 
            sym.mat[j][i] = sym.mat[i][j]; 
        }
        sym.mat[i][i] = 0.0;
    }
    return sym;
}
/* -----------------------------------------------------------------------------
 * 1.2
 * ----------------------------------------------------------------------------- */

matrix ddg(matrix data) {
    int i, j = 0;
    double sum;
    matrix matrix_sym, result_matrix;
    double **res_mat;
    if (!is_valid_matrix(&data)) {
        return empty_matrix();
    }
    matrix_sym = sym(data);
    if (!is_valid_matrix(&matrix_sym)) {
        return empty_matrix();
    }
    res_mat = init_zeros_mat(matrix_sym.rows, matrix_sym.cols);
    if(res_mat == NULL) {
        free_matrix_ptr(&matrix_sym);
        return empty_matrix();
    }
    for(i = 0; i < matrix_sym.rows; i++) {
        sum = 0;
        for (j = 0; j < matrix_sym.cols; j++) {
            sum += matrix_sym.mat[i][j];
        }
        res_mat[i][i] = sum;
    }
    result_matrix.mat = res_mat;
    result_matrix.cols = matrix_sym.cols;
    result_matrix.rows = matrix_sym.rows;
    free_matrix_ptr(&matrix_sym);
    return result_matrix;
}

/* -----------------------------------------------------------------------------
 * 1.3
 * ----------------------------------------------------------------------------- */

matrix norm(matrix data){
    matrix sym_matrix, ddg_matrix,ddg_minus_half_matrix,norm_matrix, norm_step1 ;
    sym_matrix = sym(data);
    if(!is_valid_matrix(&sym_matrix)) {
        return empty_matrix();
    }
    ddg_matrix = ddg(data);
    if(!is_valid_matrix(&ddg_matrix)) {
        free_matrixes(&sym_matrix, NULL);
        return empty_matrix();
    }
    ddg_minus_half_matrix = power_minus_half_diag(ddg_matrix);
    if(!is_valid_matrix(&ddg_minus_half_matrix)) {
        free_matrixes(&sym_matrix, &ddg_matrix, NULL);
        return empty_matrix();
    }
    norm_step1 = matrix_mult(ddg_minus_half_matrix, sym_matrix);
    if(!is_valid_matrix(&norm_step1)) {
        free_matrixes(&sym_matrix, &ddg_matrix, &ddg_minus_half_matrix, NULL);
        return empty_matrix();
    }
    norm_matrix = matrix_mult(norm_step1, ddg_minus_half_matrix);
    if(!is_valid_matrix(&norm_matrix)) {
        free_matrixes(&sym_matrix, &ddg_matrix, &ddg_minus_half_matrix,&norm_step1, NULL);
        return empty_matrix();
    }
    free_matrixes(&sym_matrix, &ddg_matrix, &ddg_minus_half_matrix, NULL);
    return norm_matrix;
}

/* -----------------------------------------------------------------------------
 * 1.4.2
 * ----------------------------------------------------------------------------- */

  matrix update_H(matrix lap_mat, matrix curr, double beta) { /*check about beta value and consider replace to const */
    int i,j;
    matrix WH,Ht,H_Ht,H_Ht_H,nextH; /* WH = W*H, Ht = H^T, H_Ht = H*H^T, H_Ht_H = H*H^T*H */
    double temp_calc = 0;
    nextH = clone_matrix(&curr);
    if(!is_valid_matrix(&nextH)) {
        return empty_matrix();
        }
    WH = matrix_mult(lap_mat,curr);
    if(!is_valid_matrix(&WH)) { 
        free_matrix_ptr(&nextH); 
        return empty_matrix(); 
        }
    Ht = transpose(curr);
    if(!is_valid_matrix(&Ht)) { 
        free_matrixes(&WH, &nextH, NULL);
        return empty_matrix(); 
        }
    H_Ht = matrix_mult(curr, Ht);
    if(!is_valid_matrix(&H_Ht)) {
        free_matrixes(&WH,&Ht,&nextH,NULL);
        return empty_matrix(); 
        }
    H_Ht_H = matrix_mult(H_Ht, curr);
    if(!is_valid_matrix(&H_Ht_H)) {
        free_matrixes(&WH,&Ht,&nextH,&H_Ht,NULL);
        return empty_matrix(); 
        }
    for(i = 0; i< curr.rows; i++) {
        for(j = 0; j< curr.cols; j++) {
            if( H_Ht_H.mat[i][j] == 0 ) { /*if denominator is zero*/
                free_matrixes(&WH,&Ht,&H_Ht_H,&H_Ht,&nextH); 
                return empty_matrix();
            }
            temp_calc = WH.mat[i][j] / H_Ht_H.mat[i][j];
            nextH.mat[i][j] = curr.mat[i][j] * (1-beta + (beta * temp_calc));
        }
    }
    free_matrixes(&WH,&Ht,&H_Ht_H,&H_Ht,NULL);
    return nextH; }


/* -----------------------------------------------------------------------------
 * 1.4.3
 * ----------------------------------------------------------------------------- */

 matrix convergence(matrix h_mat, matrix w_mat, double beta, double epsilon, int max_iter) {
    int iter = 0;
    double frobenius_norm = 0.0;
    matrix temp, nextH;
    if (!is_valid_matrix(&h_mat) || !is_valid_matrix(&w_mat)) {
        return empty_matrix();
    }

    temp = clone_matrix(&h_mat);
    if (!is_valid_matrix(&temp)) {
        return empty_matrix();
    }

    for (iter = 0; iter < max_iter; iter++) {
        nextH = update_H(w_mat, temp, beta);
        if (!is_valid_matrix(&nextH)) {
            free_matrixes(&temp, NULL);
            return empty_matrix();
        }
        frobenius_norm = squared_frobenius_norm(nextH, temp);
        if (frobenius_norm == -1) {
            free_matrixes(&temp, &nextH, NULL);
            return empty_matrix();
        }
        if (frobenius_norm < epsilon) {
            free_matrix_ptr(&temp);
            return nextH;
        }

        free_matrix_ptr(&temp);
        temp = clone_matrix(&nextH);
        free_matrix_ptr(&nextH);
        if (!is_valid_matrix(&temp)) {
            return empty_matrix();
        }
    }
    free_matrixes(&temp, &nextH, NULL);
    return empty_matrix();
}


int main(int argc, char *argv[]) {
    char *goal, *file_name;
    matrix data_mat,result_matrix;
    if(argc != 3) {
        return terminate_main(ERROR_MESSAGE);
    }
    goal = argv[1];
    file_name = argv[2];
    
    data_mat = read_data(file_name);
    if(!is_valid_matrix(&data_mat)){
        return terminate_main("read data failed");
    }
    if (strcmp(goal,"sym")==0){
        result_matrix =sym(data_mat);
    }
    else if(strcmp(goal,"ddg")==0){
        result_matrix =ddg(data_mat);
    }
    else if(strcmp(goal,"norm")==0){
        result_matrix = norm(data_mat);
    }
    else{
        free_matrixes(&data_mat, NULL);
        return terminate_main(ERROR_MESSAGE);
    }

    if(!is_valid_matrix(&result_matrix)){
        free_matrixes(&data_mat, NULL);
        return terminate_main(ERROR_MESSAGE);
    }
    print_matrix(result_matrix);
    free_matrixes(&result_matrix, &data_mat, NULL);
    return 1;
}