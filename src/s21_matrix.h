#ifndef S21_MATRIX_H_
#define S21_MATRIX_H_

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

#define S21_EPS 1e-7d
#define S21_ANTI_EPS 1e7d
#define OK 0
#define INCORRECT_MATRIX 1
#define CALC_ERROR 2

//=================create_and_remove==================//

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
void s21_reinitialization(int rows, int columns, matrix_t *target);

//======================compare========================//

#define SUCCESS 1
#define FAILURE 0

int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_is_eq_double(double, double);

//====================arithmetic=====================//

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

void s21_wide_add_matrix(int addOpt, matrix_t *, matrix_t *, matrix_t *result);

//======================transform===================//

int s21_transpose(matrix_t *A, matrix_t *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

//=====================determinant=================//
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);

void s21_minor_matrix(int, int, matrix_t *, matrix_t *);
double s21_baby_determinant(matrix_t *);
//========================math=====================//

double s21_dabs(double);
double s21_dmin(double, double);

//========================check====================//

int s21_check_correct_matrix(matrix_t *);
int s21_check_matrix_size(matrix_t *);

int s21_is_equal_size(matrix_t *, matrix_t *);

int s21_add_check(matrix_t *, matrix_t *, const matrix_t *);
int s21_mul_check(matrix_t *, matrix_t *, const matrix_t *);

int s21_check_determinant_matrix(matrix_t *);
int s21_check_inverse_matrix(matrix_t *, const matrix_t *);

#endif  // S21_MATRIX_H_
