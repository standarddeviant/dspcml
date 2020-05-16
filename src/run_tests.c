
#include <stdlib.h>
#include <stdio.h>

#define CML_IMPLEMENTATION
#include "dspcml.h"

void print_matrix(MATRIX *m, int call_number) {
	printf("Call %d:\nrows: %zd\ncols: %zd\n", call_number, m->rows, m->cols);
	for (size_t i = 0; i < m->rows; ++i) {
		printf("\t");
		for (size_t j = 0; j < m->cols; ++j) {
			printf("%6.1f", cml_get(m, i, j));
		}
		printf("\n\n");
	}
	printf("\n\n");
}

int run_cml_tests(void) {
	MATRIX *A = cml_new(3,3);
	print_matrix(A, 1);
	cml_set_all(A, 5);
	print_matrix(A, 2);
	cml_set(A, 0, 1, 99);
	print_matrix(A, 3);
	MATRIX *B = cml_identity(A->rows);
	print_matrix(B, 4);
	cml_cpy_row(B, A, 2, 0);
	print_matrix(B, 5);
	MATRIX *result = cml_new(0, 0);
	cml_mul(A, B, result);
	print_matrix(A, 6);
	print_matrix(B, 7);
	print_matrix(result, 8);
	cml_mul(A, B, NULL);
	print_matrix(A, 9);
	print_matrix(B, 10);
	

	cml_free(result);
	cml_free(B);
	cml_free(A);
	return 0;
}


int main(void) {
    int err;
    err = run_cml_tests();

    return 0;
}