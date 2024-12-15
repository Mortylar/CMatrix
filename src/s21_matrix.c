#include "s21_matrix.h"

#include <stdio.h>
#include <stdlib.h>

//====================create_and_remove=================//

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int errorStatus = OK;
  if ((rows < 1) || (columns < 1) || (!result)) {
    errorStatus = INCORRECT_MATRIX;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    for (int i = 0; i < rows; ++i) {
      result->matrix[i] = (double *)calloc(columns, sizeof(double));
    }
  }
  return errorStatus;
}

void s21_remove_matrix(matrix_t *A) {
  if (A) {
    for (int i = 0; i < A->rows; ++i) {
      free(A->matrix[i]);
      A->matrix[i] = NULL;
    }
    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
  return;
}

void s21_reinitialization(int rows, int columns, matrix_t *target) {
  s21_remove_matrix(target);
  s21_create_matrix(rows, columns, target);
  return;
}

//====================compare=====================//

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int isEqual = FAILURE;
  int errorStatus = s21_check_correct_matrix(A) || s21_check_correct_matrix(B);
  if (A && B && (A->rows == B->rows) &&
      (A->columns == B->columns) & (!errorStatus)) {
    isEqual = SUCCESS;
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; ++j) {
        if (!s21_is_eq_double(A->matrix[i][j], B->matrix[i][j])) {
          isEqual = FAILURE;
          j = A->columns;
          i = A->rows;
        }
      }
    }
  }
  return isEqual;
}

int s21_is_eq_double(double first, double second) {
  int isEqual = 0;
  if ((s21_dabs(first) < 1.0d) || (s21_dabs(second) < 1.0d)) {
    isEqual = (s21_dabs(first - second) <= S21_EPS);
  } else {
    isEqual = ((s21_dabs((first - second)) * S21_ANTI_EPS) <=
               s21_dmin(s21_dabs(first), s21_dabs(second)));
  }
  return isEqual;
}

//===================arithmetic======================//

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (result) s21_create_matrix(1, 1, result);
  int errorStatus = s21_add_check(A, B, result);
  if (!errorStatus) {
    s21_wide_add_matrix(0, A, B, result);
  }
  return errorStatus;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (result) s21_create_matrix(1, 1, result);
  int errorStatus = s21_add_check(A, B, result);
  if (!errorStatus) {
    s21_wide_add_matrix(1, A, B, result);
  }
  return errorStatus;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (result) s21_create_matrix(1, 1, result);
  int errorStatus = s21_mul_check(A, B, result);
  if (!errorStatus) {
    s21_reinitialization(A->rows, B->columns, result);
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < B->columns; ++j) {
        for (int k = 0; k < A->columns; ++k) {
          result->matrix[i][j] += (A->matrix[i][k] * B->matrix[k][j]);
        }
      }
    }
  }
  return errorStatus;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (result) s21_create_matrix(1, 1, result);
  int errorStatus = s21_check_correct_matrix(A);
  if (errorStatus || (!result)) {
    errorStatus = INCORRECT_MATRIX;
  } else {
    s21_reinitialization(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; ++j) {
        result->matrix[i][j] = number * A->matrix[i][j];
      }
    }
  }
  return errorStatus;
}

void s21_wide_add_matrix(int addOpt, matrix_t *A, matrix_t *B,
                         matrix_t *result) {
  s21_reinitialization(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; ++i) {
    for (int j = 0; j < A->columns; ++j) {
      result->matrix[i][j] = A->matrix[i][j];
      result->matrix[i][j] += ((addOpt ? (-1) : 1) * B->matrix[i][j]);
    }
  }
  return;
}

//====================transform======================//

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (result) s21_create_matrix(1, 1, result);
  int errorStatus = s21_check_correct_matrix(A);
  if ((!result) || errorStatus) {
    errorStatus = INCORRECT_MATRIX;
  } else {
    s21_reinitialization(A->columns, A->rows, result);
    for (int i = 0; i < result->rows; ++i) {
      for (int j = 0; j < result->columns; ++j) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    }
  }
  return errorStatus;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {  // A^(-1)
  if (result) s21_create_matrix(1, 1, result);
  int errorStatus = s21_check_inverse_matrix(A, result);
  if (!errorStatus) {
    s21_reinitialization(A->rows, A->columns, result);
    double detA = 0;
    errorStatus = s21_determinant(A, &detA);
    if (s21_is_eq_double(detA, 0.0)) {
      errorStatus = CALC_ERROR;
    } else {
      s21_remove_matrix(result);
      matrix_t tmpMatrix = {NULL, 0, 0};
      matrix_t antiTmpMatrix = tmpMatrix;
      s21_calc_complements(A, &tmpMatrix);
      errorStatus = (errorStatus) ? errorStatus
                                  : s21_transpose(&tmpMatrix, &antiTmpMatrix);
      errorStatus = (errorStatus) ? errorStatus
                                  : s21_mult_number(&antiTmpMatrix,
                                                    (1.0d / detA), result);
      s21_remove_matrix(&tmpMatrix);
      s21_remove_matrix(&antiTmpMatrix);
    }
  }
  return errorStatus;
}

//===================determinant=====================//

void s21_minor_matrix(int row, int column, matrix_t *daddy, matrix_t *baby) {
  s21_reinitialization(daddy->rows - 1, daddy->columns - 1, baby);
  for (int i = 0; i < daddy->rows - 1; ++i) {
    for (int j = 0; j < daddy->columns - 1; ++j)
      if ((i < row) && (j < column)) {
        baby->matrix[i][j] = daddy->matrix[i][j];
      } else if ((i < row) && (j >= column)) {
        baby->matrix[i][j] = daddy->matrix[i][j + 1];
      } else if ((i >= row) && (j < column)) {
        baby->matrix[i][j] = daddy->matrix[i + 1][j];
      } else {
        baby->matrix[i][j] = daddy->matrix[i + 1][j + 1];
      }
  }
  return;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (result) s21_create_matrix(1, 1, result);
  int errorStatus = s21_check_correct_matrix(A);
  if (!result || errorStatus) {
    errorStatus = ((result) ? errorStatus : INCORRECT_MATRIX);
  } else if (A->rows != A->columns) {
    errorStatus = CALC_ERROR;
  } else if (A->rows == 1) {
    s21_reinitialization(1, 1, result);
    result->matrix[0][0] = 1.0d;
  } else {
    int mSize = A->rows;
    s21_reinitialization(mSize, mSize, result);
    double localDet = 0;
    matrix_t localMatrix = {NULL, 0, 0};
    for (int i = 0; i < mSize; ++i) {
      for (int j = 0; j < mSize; ++j) {
        s21_minor_matrix(i, j, A, &localMatrix);
        s21_determinant(&localMatrix, &localDet);
        s21_remove_matrix(&localMatrix);
        int doubtingUnit = ((((i + j) % 2)) ? (-1) : 1);
        result->matrix[i][j] = (doubtingUnit * localDet);
      }
    }
    s21_remove_matrix(&localMatrix);
  }
  return errorStatus;
}

int s21_determinant(matrix_t *A, double *result) {
  int errorStatus = s21_check_determinant_matrix(A);
  errorStatus = (((!result) && (!errorStatus)) ? CALC_ERROR : errorStatus);
  if ((!errorStatus) && (result)) {
    *result = 0;
    if (A->rows > 3) {
      int doubtingUnit = 1;
      matrix_t thisMinorMatrix = {NULL, 0, 0};
      for (int i = 0; i < A->rows; ++i) {
        double tmpMinor = 0;
        s21_minor_matrix(i, 0, A, &thisMinorMatrix);
        s21_determinant(&thisMinorMatrix, &tmpMinor);
        *result += (doubtingUnit * tmpMinor * A->matrix[i][0]);
        doubtingUnit *= (-1);
      }
      s21_remove_matrix(&thisMinorMatrix);
    } else {
      *result = s21_baby_determinant(A);
    }
  }
  return errorStatus;
}

double s21_baby_determinant(matrix_t *A) {
  double result = 0;
  int mSize = A->rows;
  if (1 == mSize) {
    result = A->matrix[0][0];
  } else if (2 == mSize) {
    result = (A->matrix[0][0] * A->matrix[1][1]);
    result -= (A->matrix[1][0] * A->matrix[0][1]);
  } else if (3 == mSize) {
    result = (A->matrix[0][0] * A->matrix[1][1] * A->matrix[2][2]);
    result += (A->matrix[0][1] * A->matrix[1][2] * A->matrix[2][0]);
    result += (A->matrix[0][2] * A->matrix[1][0] * A->matrix[2][1]);
    result -= (A->matrix[2][0] * A->matrix[1][1] * A->matrix[0][2]);
    result -= (A->matrix[0][0] * A->matrix[2][1] * A->matrix[1][2]);
    result -= (A->matrix[1][0] * A->matrix[0][1] * A->matrix[2][2]);
  }
  return result;
}

//======================math==========================//

double s21_dabs(double number) { return (number < 0) ? (-number) : number; }

double s21_dmin(double first, double second) {
  return ((first < second) ? first : second);
}

//======================check=========================//

int s21_check_correct_matrix(matrix_t *matrix) {
  int errorStatus = INCORRECT_MATRIX;
  if (matrix) {
    errorStatus = s21_check_matrix_size(matrix);
  }
  return errorStatus;
}

int s21_check_matrix_size(matrix_t *matrix) {
  int errorStatus = OK;
  if ((matrix->rows < 1) || (matrix->columns < 1)) {
    errorStatus = INCORRECT_MATRIX;
  }
  return errorStatus;
}

int s21_is_equal_size(matrix_t *first, matrix_t *second) {
  return ((first->rows == second->rows) && (first->columns == second->columns));
}

int s21_add_check(matrix_t *A, matrix_t *B, const matrix_t *result) {
  int errorStatus =
      (s21_check_correct_matrix(A) || s21_check_correct_matrix(B));
  if (!result) {
    errorStatus = INCORRECT_MATRIX;
  } else if (!s21_is_equal_size(A, B) && (!errorStatus)) {
    errorStatus = CALC_ERROR;
  }
  return errorStatus;
}

int s21_mul_check(matrix_t *A, matrix_t *B, const matrix_t *result) {
  int errorStatus =
      (s21_check_correct_matrix(A) || s21_check_correct_matrix(B));
  if (!result) {
    errorStatus = INCORRECT_MATRIX;
  } else if ((A->columns != B->rows) && (!errorStatus)) {
    errorStatus = CALC_ERROR;
  }
  return errorStatus;
}

int s21_check_determinant_matrix(matrix_t *A) {
  int errorStatus = s21_check_correct_matrix(A);
  if (!errorStatus) {
    errorStatus = ((A->rows == A->columns) ? errorStatus : CALC_ERROR);
  }
  return errorStatus;
}

int s21_check_inverse_matrix(matrix_t *A, const matrix_t *result) {
  int errorStatus = s21_check_determinant_matrix(A);
  if (!errorStatus) {
    if (!result) {
      errorStatus = INCORRECT_MATRIX;
    }
  }
  return errorStatus;
}
