
#include <stdio.h>
#include <errno.h>

#define CML_NO_DEPENDENCIES
#define CML_IMPLEMENTATION
#include "dspcml.h"



// #define STRERR_BUF_LEN (128)
// char g_strerr_buf[STRERR_BUF_LEN];
size_t maybe_print_errno(char *cptr) {
    cptr[0] = '\0';
    if (errno) {
        strerror_r(errno, cptr, 4096);
    }
    return (size_t)strlen(cptr);
}


size_t print_matrix(unsigned long ulptr, char *cptr) {
    size_t length = 0;
    MATRIX *m = (MATRIX *)ulptr;
	length += sprintf(cptr + length,
        "MATRIX: sz = %zd x  %zd\n",
        m->rows, m->cols);
	for (size_t i = 0; i < m->rows; ++i) {
		length += sprintf(cptr + length, "    ");
		for (size_t j = 0; j < m->cols; ++j) {
			length += sprintf(cptr + length, "%e  ", cml_get(m, i, j));
		}
		length += sprintf(cptr + length, "\n");
	}
	length += sprintf(cptr + length, "\n");

    return length;
}


unsigned long zeros(size_t rows, size_t cols) {
    return (unsigned long)cml_zeros(rows, cols);
}


unsigned long arange(cml_real_t start, cml_real_t stop, cml_real_t step) {
    return (unsigned long)cml_arange(start, stop, step);
}


unsigned long tile(unsigned long ulptr, size_t row_copies, size_t col_copies) {
    return (unsigned long)cml_tile((MATRIX *)ulptr, row_copies, col_copies);
}


unsigned long new_sos(size_t nb_bq, size_t nb_ch, cml_real_t *coeffs) {
    return (unsigned long )cml_new_sos(nb_bq, nb_ch, coeffs);
}


void cpy_submat(unsigned long ulptr_dst, unsigned long ulptr_src, size_t dst_row, size_t dst_col, size_t src_row, size_t src_col, size_t nb_rows, size_t nb_cols) {
    cml_cpy_submat(
        (MATRIX *)ulptr_dst,
        (MATRIX *)ulptr_src,
        dst_row, dst_col,
        src_row, src_col,
        nb_rows, nb_cols
    );
}


void __sin__(unsigned long ulptr, unsigned long ulptr_opt) {
    cml_sin((MATRIX *)ulptr, (MATRIX *)ulptr_opt);
}


void __cos__(unsigned long ulptr, unsigned long ulptr_opt) {
    cml_cos((MATRIX *)ulptr, (MATRIX *)ulptr_opt);
}

void reshape(unsigned long ulptr, size_t new_rows, size_t new_cols) {
    cml_reshape((MATRIX *)ulptr, new_rows, new_cols);
}


size_t numel(unsigned long ulptr) {
    MATRIX *m = (MATRIX *)ulptr;
    /* TODO checks on m */

    return (unsigned long)(m->rows * m->cols);
}


size_t rows(unsigned long ulptr) {
    MATRIX *m = (MATRIX *)ulptr;
    /* TODO checks on m */

    return (unsigned long)m->rows;
}


size_t cols(unsigned long ulptr) {
    MATRIX *m = (MATRIX *)ulptr;
    /* TODO checks on m */

    return (unsigned long)m->cols;
}


void to_numpy(unsigned long ulptr, cml_real_t *out_data) {
    MATRIX *m = (MATRIX *)ulptr;
    /* TODO checks on m */

    printf("rows = %lu\n", m->rows);
    printf("cols = %lu\n", m->cols);
    
    /* copy cml_real_t data from matrix to out_data */
    memcpy(
        (void *)out_data,                       /* dst      */
        (const void *)m->data,                  /* src      */
        m->rows * m->cols * sizeof(cml_real_t)  /* nb bytes */
    );
}


unsigned long to_dspcml(cml_real_t *in_data, size_t rows, size_t cols) {
    /* TODO checks... */
    MATRIX *m = cml_new(rows, cols);

    memcpy( 
        (void *)m->data,                        /* src      */
        (const void *)in_data,                  /* dst      */
        m->rows * m->cols * sizeof(cml_real_t)  /* nb bytes */
    );

    return (unsigned long)m;
}



