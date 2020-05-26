
#include <stdio.h>
#include <errno.h>

#define CML_NO_DEPENDENCIES
#define CML_IMPLEMENTATION
#include "dspcml.h"



#define STRERR_BUF_LEN (128)
char g_strerr_buf[STRERR_BUF_LEN];
void maybe_print_errno() {
    if (errno) {
        strerror_r(errno, g_strerr_buf, STRERR_BUF_LEN);
        fprintf(stdout, "errno = %s\n", g_strerr_buf);
    }
}


unsigned int print_matrix(unsigned long ulptr, char *cptr) {
    size_t length = 0;
    MATRIX *m = (MATRIX *)ulptr;
	length += sprintf(cptr + length,
        "MATRIX: sz = %zd x  %zd\n",
        m->rows, m->cols);
	for (size_t i = 0; i < m->rows; ++i) {
		length += sprintf(cptr + length, "    ");
		for (size_t j = 0; j < m->cols; ++j) {
			length += sprintf(cptr + length, "%6.1f", cml_get(m, i, j));
		}
		length += sprintf(cptr + length, "\n");
	}
	length += sprintf(cptr + length, "\n\n");

    return length;
}


unsigned long zeros(size_t rows, size_t cols) {
    void *ptr = (void *)cml_zeros(rows, cols);
    maybe_print_errno();
    return (unsigned long)ptr;
}


unsigned long numel(unsigned long ulptr) {
    MATRIX *m = (MATRIX *)ulptr;
    /* TODO checks on m */

    return (unsigned long)(m->rows * m->cols);
}


unsigned long rows(unsigned long ulptr) {
    MATRIX *m = (MATRIX *)ulptr;
    /* TODO checks on m */

    return (unsigned long)m->rows;
}


unsigned long cols(unsigned long ulptr) {
    MATRIX *m = (MATRIX *)ulptr;
    /* TODO checks on m */

    return (unsigned long)m->cols;
}


void to_numpy(unsigned long ulptr, float *out_data) {
    MATRIX *m = (MATRIX *)ulptr;
    /* TODO checks on m */

    printf("rows = %lu\n", m->rows);
    printf("cols = %lu\n", m->cols);
    
    /* copy cml_real_t data from matrix to out_data */
    memcpy(
        (void *)out_data,                       /* dst      */
        (const void *)m->data,                  /* src      */
        m->rows * m->cols * sizeof(float)  /* nb bytes */
    );
}




