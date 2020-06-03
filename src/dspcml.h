/* CML 1.0.0 -- A header-only C matrix library.
 *
 * ABOUT: 
 *      CML is a (fully C++ compatible) header-only C matrix library designed with a focus on portability, simplicity, 
 *      and efficiency written under the POSIX standard. Several computationally-intense functions require/are enhanced 
 *      by BLAS/LAPACK. The user is able to define optional flags prior to including this library in source code to do 
 *      things such as removing BLAS/LAPACK as a dependency (while simultaneously removing several functions from API), 
 *   	or giving all functions static storage. CML operates on type 'MATRIX' for all processes. 'MATRIX' is column major 
 *   	so the user can easily interface with the most popular linear-algebra-related API's. 'errno' is used for error 
 *   	handling extensively. 'errno' is set only in the presence of an error. The user is free to disregard 'errno' if 
 *   	an error occurred, but the state of the resultant matrix is undefined. 
 *
 * USAGE:
 *      CML can be used in either header-only mode or implementation mode. The header-only mode is used by default when
 *      included and allows including this header in other headers and does not contain the actual implementation. The
 *      implementation mode requires the user to define the preprocessor macro 'CML_IMPLEMENTATION' in one .c/.cpp file
 *      before #including this file, e.g.:
 *          
 *          #define CML_IMPLEMENTATION
 *          #include "cml.h"
 *
 *      Also optionally define the flags listed in the section 'OPTIONAL DEFINES' below in header and implementation
 *      mode if you want to use additional functionality or need more control over the library. Everytime you include 
 *      "cml.h", you have to define the same optional flags. Not doing so will lead to compile-time errors.
 *
 * OPTIONAL DEFINES:
 *      CML_PRIVATE
 *          If defined, declares all functions as static, so they can only be accessed inside the file that contains the
 *          implementation.
 *
 *      CML_NO_DEPENDENCIES
 *          If defined, will remove all functions that fully depend on external API's. The user will then be able to
 *          compile and link without the external API's.
 */
#ifndef CML_H
#define CML_H


#ifdef __cplusplus
extern "C" {
#endif


#include <stddef.h>
#include <stdbool.h>


#ifdef CML_PRIVATE
    #define CML_API static inline
#else
    #define CML_API extern
#endif


#include "dspcml_element_type.h"

/* ========================================================
 *
 *                          API
 *
 * ======================================================== */
typedef struct {
    size_t rows, cols;
    cml_real_t *data;
} MATRIX;
typedef struct ZMATRIX {
    size_t rows, cols;
    cml_cpx_t *data;
} ZMATRIX;
#define MATRIX_COL2PTR(m, cix) ( (cml_real_t *)&(m->data[cix*m->rows]) )
#define MATRIX_CHAN2PTR(m, cix) ( (cml_real_t *)&(m->data[cix*m->rows]) )
#define MATRIX_ALL_CH (SIZE_MAX) /* like the colon operator in NumPy */
/* stride length to iterate across a certain dim */
#define MATRIX_STRIDE_DIM(m, d) ( 0==d ? 1 : m->rows )


/*    MATRIX SIZE CHECK MACROS    */
#define DIMSZ(m, d) (0 == d ? m->rows : m->cols)
#define SAME_DIM_SZ(m1, d1, m2, d2) (DIMSZ(m1,d1) == DIMSZ(m2,d2))
#define DIFF_DIM_SZ(m1, d1, m2, d2) (DIMSZ(m1,d1) != DIMSZ(m2,d2))
#define SAME_DIMS(m1, m2) ((DIMSZ(m1, 0) == DIMSZ(m2, 0)) && (DIMSZ(m1, 1) == DIMSZ(m2, 1)))
#define DIFF_DIMS(m1, m2) ((DIMSZ(m1, 0) != DIMSZ(m2, 0)) || (DIMSZ(m1, 1) != DIMSZ(m2, 1))) 


/*    ALLOCATION    */
CML_API  MATRIX*       cml_new(size_t rows, size_t cols);
CML_API  MATRIX*       cml_dup(MATRIX *m);
CML_API  MATRIX*       cml_zeros(size_t rows, size_t cols);
CML_API  MATRIX*       cml_ones(size_t rows, size_t cols);
CML_API  MATRIX*       cml_zeros_like(MATRIX *m);
CML_API  MATRIX*       cml_ones_like(MATRIX *m);
CML_API  MATRIX*       cml_tile(MATRIX *m, size_t row_copies, size_t col_copies);
CML_API  MATRIX*       cml_identity(size_t dim);
CML_API  MATRIX*       cml_lower_tri(size_t dim);
CML_API  MATRIX*       cml_upper_tri(size_t dim);
CML_API  MATRIX*       cml_arange(cml_real_t start, cml_real_t stop, cml_real_t step);
CML_API  MATRIX*       cml_new_sos(size_t nb_bq, size_t nb_ch, cml_real_t *coeffs);
CML_API  MATRIX*       cml_new_fir_mat(size_t taps, size_t len, size_t chans);
CML_API  ZMATRIX*      cml_real_to_cpx(MATRIX *m);
CML_API  void          cml_free(MATRIX *m);


/*    GET/SET    */
CML_API  cml_real_t    cml_get(MATRIX *m, size_t row, size_t col);
CML_API  void          cml_set(MATRIX *m, size_t row, size_t col, cml_real_t data);
CML_API  void          cml_set_all(MATRIX *m, cml_real_t data);
CML_API  void          cml_set_row(MATRIX *m, size_t row, cml_real_t data);
CML_API  void          cml_set_col(MATRIX *m, size_t col, cml_real_t data);
CML_API  void          cml_set_sos_coeffs_all(MATRIX *sos, cml_real_t *coeffs);
CML_API  void          cml_reshape(MATRIX *m, size_t new_rows, size_t new_cols);


/*    PHYSICAL MANIPULATION    */
CML_API  void          cml_cpy(MATRIX *dst, MATRIX *src);
CML_API  void          cml_cpy_row(MATRIX *dst, MATRIX *src, size_t dst_row, size_t src_row);
CML_API  void          cml_cpy_row_range(MATRIX *dst, MATRIX *src, size_t dst_row, size_t src_row, size_t nb_rows);
CML_API  void          cml_cpy_col(MATRIX *dst, MATRIX *src, size_t dst_col, size_t src_col);
CML_API  void          cml_cpy_col_range(MATRIX *dst, MATRIX *src, size_t dst_col, size_t src_col, size_t nb_cols);
CML_API  void          cml_cpy_elem(MATRIX *dst, MATRIX *src, size_t dst_row, size_t dst_col, size_t src_row, size_t src_col);
CML_API  void          cml_cpy_submat(MATRIX *dst, MATRIX *src, size_t dst_row, size_t dst_col, size_t src_row, size_t src_col, size_t nb_rows, size_t nb_cols); 
CML_API  void          cml_cpy_self_row(MATRIX *m, size_t dst_row, size_t src_row);
CML_API  void          cml_cpy_self_col(MATRIX *m, size_t dst_col, size_t src_col);
CML_API  void          cml_cpy_self_elem(MATRIX *m, size_t dst_row, size_t dst_col, size_t src_row, size_t src_col);
CML_API  void          cml_cpy_self_row_range(MATRIX *m, size_t dst_row, size_t src_row, size_t nb_rows);
CML_API  void          cml_cpy_self_col_range(MATRIX *m, size_t dst_col, size_t src_col, size_t nb_cols);
CML_API  void          cml_swap(MATRIX *m1, MATRIX *m2);
CML_API  void          cml_swap_row(MATRIX *m1, MATRIX *m2, size_t m1_row, size_t m2_row);
CML_API  void          cml_swap_col(MATRIX *m1, MATRIX *m2, size_t m1_col, size_t m2_col);
CML_API  void          cml_swap_elem(MATRIX *m1, MATRIX *m2, size_t m1_row, size_t m1_col, size_t m2_row, size_t m2_col);
CML_API  void          cml_swap_self_row(MATRIX *m, size_t row_1, size_t row_2);
CML_API  void          cml_swap_self_col(MATRIX *m, size_t col_1, size_t col_2);
CML_API  void          cml_swap_self_elem(MATRIX *m, size_t row_1, size_t col_1, size_t row_2, size_t col_2);
CML_API  void          cml_ins_row(MATRIX *dst, MATRIX *src, size_t dst_row, size_t src_row);
CML_API  void          cml_ins_col(MATRIX *dst, MATRIX *src, size_t dst_col, size_t src_col);
CML_API  void          cml_ins_self_row(MATRIX *m, size_t dst_row, size_t src_row);
CML_API  void          cml_ins_self_col(MATRIX *m, size_t dst_col, size_t src_col);
CML_API  void          cml_del_all(MATRIX *m);
CML_API  void          cml_del_row(MATRIX *m, size_t row);
CML_API  void          cml_del_col(MATRIX *m, size_t col);
CML_API  void          cml_adjoin_top(MATRIX *dst, MATRIX *src);
CML_API  void          cml_adjoin_bottom(MATRIX *dst, MATRIX *src);
CML_API  void          cml_adjoin_left(MATRIX *dst, MATRIX *src);
CML_API  void          cml_adjoin_right(MATRIX *dst, MATRIX *src);


/*    ARITHMETIC    */
CML_API  void          cml_add_const(MATRIX *m, cml_real_t data, MATRIX *opt);
CML_API  void          cml_mul_const(MATRIX *m, cml_real_t data, MATRIX *opt);
CML_API  void          cml_add(MATRIX *m1, MATRIX *m2, MATRIX *opt);
CML_API  void          cml_sub(MATRIX *m1, MATRIX *m2, MATRIX *opt);
CML_API  void          cml_mul(MATRIX *m1, MATRIX *m2, MATRIX *opt);
CML_API  void          cml_mul_elem(MATRIX *m1, MATRIX *m2, MATRIX *opt);
CML_API  void          cml_div_elem(MATRIX *m1, MATRIX *m2, MATRIX *opt);


/* TRIGONOMETRY */
CML_API void cml_sin(MATRIX *m, MATRIX *opt);
CML_API void cml_cos(MATRIX *m, MATRIX *opt);



/*    SIGNAL PROCESSSING (MATRIX ONLY)   */
CML_API  void          cml_add_chan_m1_m2(MATRIX *m1, size_t m1ch, MATRIX *m2, size_t m2ch);
CML_API  cml_real_t    cml_sum(MATRIX *m);
CML_API  void          cml_sum_dim(MATRIX *m, size_t dim, MATRIX *opt);
CML_API  cml_real_t    cml_mean(MATRIX *m);
CML_API  void          cml_mean_dim(MATRIX *m, size_t dim, MATRIX *opt);
CML_API  cml_real_t    cml_rms(MATRIX *m);
CML_API  void          cml_rms_dim(MATRIX *m, size_t dim, MATRIX *opt);
CML_API  void          cml_abs(MATRIX *m, MATRIX *opt);
CML_API  void          cml_pow_mat2const(MATRIX *m, cml_real_t r1, MATRIX *opt);
CML_API  void          cml_pow_const2mat(cml_real_t r1, MATRIX *m, MATRIX *opt);
CML_API  void          cml_10_log10_abs(MATRIX *m, MATRIX *opt);
CML_API  void          cml_20_log10_abs(MATRIX *m, MATRIX *opt);
CML_API  void          cml_clip(MATRIX *m, cml_real_t lolim, cml_real_t hilim, MATRIX *opt);
CML_API  void          cml_sos_proc(MATRIX *sos, MATRIX *m, MATRIX *scratch, MATRIX *opt);



/*    OTHER OPERATIONS    */
CML_API  cml_real_t    cml_min(MATRIX *m);
CML_API  cml_real_t    cml_max(MATRIX *m);
CML_API  void          cml_min_dim(MATRIX *m, size_t dim, MATRIX *opt);
CML_API  void          cml_max_dim(MATRIX *m, size_t dim, MATRIX *opt);
CML_API  bool          cml_is_zero(MATRIX *m);
CML_API  bool          cml_is_pos(MATRIX *m);
CML_API  bool          cml_is_neg(MATRIX *m);
CML_API  bool          cml_is_nonneg(MATRIX *m);
CML_API  bool          cml_is_equal(MATRIX *m1, MATRIX *m2);
CML_API  void          cml_transpose(MATRIX *m, MATRIX *opt);
CML_API  void          cml_normalize(MATRIX *m, MATRIX *opt);


/*    OTHER OPERATIONS (BLAS/LAPACK REQUIRED)    */
#ifndef CML_NO_DEPENDENCIES
CML_API  bool          cml_inverse(MATRIX *m, MATRIX *opt);
CML_API  bool          cml_sys_equ(MATRIX *A, MATRIX *B, MATRIX *opt);
#endif


#ifdef __cplusplus
}
#endif


/* ========================================================
 *
 *                      IMPLEMENTATION
 *
 * ======================================================== */
#ifdef CML_IMPLEMENTATION


#ifdef __cplusplus
extern "C" {
#endif


#ifdef CML_PRIVATE
    #define VEC_PRIVATE
#endif


#include <stdlib.h>
#include <errno.h>
#include "vec.h" 


#ifndef CML_NO_DEPENDENCIES
/* TODO - handle single vs. double routines */ 
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
void dgetrf_(int*, int*, double*, int*, int*, int*);
void dgetri_(int*, double*, int*, int*, double*, int*, int*);
void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
#endif


#ifdef __cplusplus
}
#endif


CML_API MATRIX* cml_new(size_t rows, size_t cols) {
    if ((rows == 0 && cols > 0) || (cols == 0 && rows > 0)) {
        errno = EINVAL;
        return NULL;
    }

    int old_errno = errno;
    errno = 0;
    MATRIX *m = (MATRIX *) malloc(sizeof(MATRIX));

    if (m == NULL) {
        return NULL;
    }

    m->rows = rows;
    m->cols = cols;
    m->data = (cml_real_t *) vec_new_cap(sizeof(cml_real_t), rows * cols);

    if (m->data == NULL) {
        free(m);
        return NULL;
    } else {
        errno = old_errno;
    }

    for (size_t i = 0; i < rows * cols; ++i) {
        vec_push(m->data, 0);
    }

    return m;
}


CML_API MATRIX* cml_dup(MATRIX *m) {
    if (m == NULL) {
        errno = EINVAL;
        return NULL;
    }

    MATRIX *dup = cml_new(m->rows, m->cols);

    if (dup == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < dup->rows; ++i) {
        for (size_t j = 0; j < dup->cols; ++j) {
            cml_set(dup, i, j, cml_get(m, i, j));
        }
    }

    return dup;
}

CML_API MATRIX* cml_zeros(size_t rows, size_t cols) {
    return cml_new(rows, cols);
}


CML_API MATRIX* cml_ones(size_t rows, size_t cols) {
    MATRIX *m = cml_new(rows, cols);

    if (m == NULL) {
        return NULL;
    }

    cml_set_all(m, 1);

    return m;
}


CML_API  MATRIX* cml_zeros_like(MATRIX *m) {
    return cml_new(m->rows, m->cols);
}


CML_API MATRIX* cml_ones_like(MATRIX *m) {
    return cml_ones(m->rows, m->cols);
}


CML_API MATRIX* cml_tile(MATRIX *m, size_t row_copies, size_t col_copies) {
    /* TODO checks... */
    MATRIX *opt = cml_new(m->rows * row_copies, m->cols * col_copies);

    for(size_t rc = 0; rc < row_copies; ++rc) {
        for(size_t cc = 0; cc < col_copies; ++cc) {
            cml_cpy_submat(opt, m, rc * m->rows, cc * m->cols, 0, 0, m->rows, m->cols);
        }
    }

    return opt;
}


CML_API MATRIX* cml_identity(size_t dim) {
    MATRIX *m = cml_new(dim, dim);

    if (m == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < dim; ++i) {
        cml_set(m, i, i, 1);
    }

    return m;
}


CML_API MATRIX* cml_lower_tri(size_t dim) {
    MATRIX *m = cml_new(dim, dim);

    if (m == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            cml_set(m, i, j, 1);
        }
    }

    return m;
}


CML_API MATRIX* cml_upper_tri(size_t dim) {
    MATRIX *m = cml_new(dim, dim);

    if (m == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = i; j < dim; ++j) {
            cml_set(m, i, j, 1);
        }
    }

    return m;
}


/* execute scalar function on each element with no arguments per element */
CML_API void cml_elm_real_func0(MATRIX *m, cml_real_func0_t op, MATRIX *opt) {
    /* TODO - checks on m / opt */
    if(NULL == opt) {
        opt = m;
    }
    for (size_t i = 0; i < vec_len(m->data); ++i) {
        opt->data[i] = op(m->data[i]);
    }
}


/* execute scalar function on each element with 1 argument per element */
CML_API void cml_elm_real_func1(MATRIX *m, cml_real_func1_t op, cml_real_t r1, MATRIX *opt) {
    /* TODO - checks on m / opt */
    if(NULL == opt) {
        opt = m;
    }
    for (size_t i = 0; i < vec_len(m->data); ++i) {
        opt->data[i] = op(m->data[i], r1);
    }
}


/* reverse call signature order btween m and r1 compared to cml_elm_real_func1 */
CML_API void cml_elm_real_1func(cml_real_t r1, cml_real_func1_t op, MATRIX *m, MATRIX *opt) {
    /* TODO - checks on m / opt */
    if(NULL == opt) {
        opt = m;
    }
    for (size_t i = 0; i < vec_len(m->data); ++i) {
        opt->data[i] = op(r1, m->data[i]);
    }
}


/* execute scalar function on each element with 2 argument per element */
CML_API void cml_elm_real_func2(MATRIX *m, cml_real_func2_t op, cml_real_t r1, cml_real_t r2, MATRIX *opt) {
    /* TODO - checks on m / opt */
    if(NULL == opt) {
        opt = m;
    }
    for (size_t i = 0; i < vec_len(m->data); ++i) {
        opt->data[i] = op(m->data[i], r1, r2);
    }
}


CML_API void cml_cos(MATRIX *m, MATRIX *opt) {
    /* TODO - checks on m / opt */
    cml_elm_real_func0(m, cml_real_cos, opt);
}


CML_API void cml_sin(MATRIX *m, MATRIX *opt) {
    /* TODO - checks on m / opt */
    cml_elm_real_func0(m, cml_real_sin, opt);
}


/* linear generator */
CML_API MATRIX* cml_arange(cml_real_t start, cml_real_t stop, cml_real_t step) {
    cml_real_t tmp = (stop - start) / step;
    size_t rows = (size_t)cml_real_ceil(tmp);

    MATRIX *m = cml_new(rows, 1);
    if(NULL == m) {
        return m;
    }

    tmp = start;
    for (size_t i = 0; i < rows; ++i) {
        cml_set(m, i, 0, tmp);
        tmp += step;
    }

    return m;
}


// CML_API MATRIX* cml0_gen_cos(size_t rows, size_t cols, cml_real_t phase_incr, MATRIX *phase_starts) {
//     /* create matrix of phases values */
//     MATRIX *m = cml0_gen_lin(rows, cols, phase_incr, phase_starts);
//     cml_cos(m, NULL); /* cos of phases in-place */
//     return m;
// }


// CML_API MATRIX* cml0_gen_sin(size_t rows, size_t cols, cml_real_t phase_incr, MATRIX *phase_starts) {
//     /* create matrix of phases values */
//     MATRIX *m = cml0_gen_lin(rows, cols, phase_incr, phase_starts);
//     cml_sin(m, NULL); /* sin of phases in-place */
//     return m;
// }


/* multi-channel SOS filter format that stores the coefficients with the state 
   'back-to-back' and per-channel. This leverages the cml MATRIX to store
   multi-channel SOS or 'cascaded biqauds' data in one place and support
   multi-channel filtering with no extra data types. A matrix of this format
   could be constructed directly and manually */
CML_API MATRIX* cml_new_sos(size_t nb_bq, size_t nb_ch, cml_real_t *coeffs) {
    MATRIX *sos = cml_new(7*nb_bq, nb_ch);
    if (sos == NULL) {
        return NULL;
    }

    for (size_t ch = 0; ch < nb_ch; ++ch) {
        for (size_t bq = 0; bq < nb_bq; ++bq) {
            cml_real_t *bq_coef = coeffs + 6 * (bq + (ch * nb_bq));
            cml_set(sos, 7*bq + 0, ch, *(bq_coef + 0)); /* b0 */
            cml_set(sos, 7*bq + 1, ch, *(bq_coef + 1)); /* b1 */
            cml_set(sos, 7*bq + 2, ch, *(bq_coef + 2)); /* b2 */
            cml_set(sos, 7*bq + 3, ch, *(bq_coef + 4)); /* a1 */
            cml_set(sos, 7*bq + 4, ch, *(bq_coef + 5)); /* a2 */
            cml_set(sos, 7*bq + 5, ch, 0);              /* w0 */
            cml_set(sos, 7*bq + 6, ch, 0);              /* w1 */
        }
    }

    return sos;
}


/* very simple function to create new fir filter state */
CML_API MATRIX* cml_new_fir_state(size_t taps, size_t len, size_t chans) {
    return cml_new(taps + len - 1, chans);
}


CML_API ZMATRIX* cml_convert_MATRIX_to_ZMATRIX(MATRIX *m) {
    ZMATRIX *z = (ZMATRIX *)m;
    z->cols = z->rows >> 1; /* this should do the trick */
    return z;
}


CML_API void cml_free(MATRIX *m) {
    if (m == NULL) {
        return;
    }

    vec_free(m->data);
    free(m);
}


/* TODO - test cmlz_free */
CML_API void cmlz_free(ZMATRIX *z) {
    cml_free((MATRIX*)z);
}


CML_API cml_real_t cml_get(MATRIX *m, size_t row, size_t col) {
    if (m == NULL || row >= m->rows || col >= m->cols) {
        errno = EINVAL;
        return 0;
    }

    return *(m->data + m->rows * col + row);
}


CML_API void cml_set(MATRIX *m, size_t row, size_t col, cml_real_t data) {
    if (m == NULL || row >= m->rows || col >= m->cols) {
        errno = EINVAL;
        return;
    }

    *(m->data + m->rows * col + row) = data;
}


CML_API void cml_set_all(MATRIX *m, cml_real_t data) {
    if (m == NULL) {
        errno = EINVAL;
        return;
    }

    for (size_t i = 0; i < vec_len(m->data); ++i) {
        m->data[i] = data;
    }
}


CML_API void cml_set_row(MATRIX *m, size_t row, cml_real_t data) {
    if (m == NULL || row >= m->rows) {
        errno = EINVAL;
        return;
    }

    for (size_t j = 0; j < m->cols; ++j) {
        cml_set(m, row, j, data);
    }
}


CML_API void cml_set_col(MATRIX *m, size_t col, cml_real_t data) {
    if (m == NULL || col >= m->cols) {
        errno = EINVAL;
        return;
    }

    for (size_t i = 0; i < m->rows; ++i) {
        cml_set(m, i, col, data);
    }
}


CML_API void cml_set_sos_state_zero(MATRIX *sos) {
    for (size_t ch = 0; ch < sos->cols; ++ch) {
        size_t row = 0;
        while(row < sos->rows) {
            row += 5;
            cml_set(sos, row++, ch, 0); /* w0 */
            cml_set(sos, row++, ch, 0); /* w1 */
        }
    }
}


/* NOTE: this does 'column-major' reshaping and does not modify m->data */
CML_API void cml_reshape(MATRIX *m, size_t new_rows, size_t new_cols) {
    if (m->rows * m->cols != new_rows * new_cols) {
        errno = EINVAL;
        return;
    }

    m->rows = new_rows;
    m->cols = new_cols;
}


CML_API void cml_cpy(MATRIX *dst, MATRIX *src) {
    if (dst == NULL || src == NULL || dst == src) {
        errno = EINVAL;
        return;
    }

    int old_errno = errno;
    errno = 0;
    cml_real_t *copied_data = (cml_real_t *) vec_dup(src->data);

    if (errno) {
        return;
    } else {
        errno = old_errno;
    }

    vec_free(dst->data);
    dst->rows = src->rows;
    dst->cols = src->cols;
    dst->data = copied_data;
}


CML_API void cml_cpy_row(MATRIX *dst, MATRIX *src, size_t dst_row, size_t src_row) {
    if (dst == NULL || src == NULL || dst == src || dst->cols != src->cols || dst_row >= dst->rows || src_row >= src->rows) {
        errno = EINVAL;
        return;
    }

    for (size_t j = 0; j < dst->cols; ++j) {
        cml_set(dst, dst_row, j, cml_get(src, src_row, j));
    }
}

CML_API void cml_cpy_row_range(MATRIX *dst, MATRIX *src, size_t dst_row, size_t src_row, size_t nb_rows) {
    if (dst == NULL || src == NULL || dst->cols != src->cols ||  dst_row + nb_rows > dst->rows || src_row + nb_rows > src->rows) {
        errno = EINVAL;
        return;
    }
    
    for(size_t j = 0; j < src->cols; ++j) {
        // void * memcpy ( void * dstination, const void * source, size_t num );
        memcpy(
            (void *)(MATRIX_COL2PTR(dst, j) + dst_row),
            (void *)(MATRIX_COL2PTR(src, j) + src_row),
            (nb_rows * sizeof(cml_real_t))
        );
    }
}



CML_API void cml_cpy_col(MATRIX *dst, MATRIX *src, size_t dst_col, size_t src_col) {
    if (dst == NULL || src == NULL || dst == src || dst->rows != src->rows || dst_col >= dst->cols || src_col >= src->cols) {
        errno = EINVAL;
        return;
    }

    for (size_t i = 0; i < dst->rows; ++i) {
        cml_set(dst, i, dst_col, cml_get(src, i, src_col));
    }
}


CML_API void cml_cpy_elem(MATRIX *dst, MATRIX *src, size_t dst_row, size_t dst_col, size_t src_row, size_t src_col) {
    if (dst == NULL || src == NULL || dst == src || dst_row >= dst->rows || dst_col >= dst->cols || src_row >= src->rows || src_col >= src->cols) {
        errno = EINVAL;
        return;
    }

    cml_set(dst, dst_row, dst_col, cml_get(src, src_row, src_col));
}


CML_API void cml_cpy_self_row(MATRIX *m, size_t dst_row, size_t src_row) {
    if (m == NULL || dst_row >= m->rows || src_row >= m->rows) {
        errno = EINVAL;
        return;
    }

    for (size_t j = 0; j < m->cols; ++j) {
        cml_set(m, dst_row, j, cml_get(m, src_row, j));
    }
}


CML_API void cml_cpy_self_row_range(MATRIX *m, size_t dst_row, size_t src_row, size_t nb_rows) {
    cml_cpy_row_range(m, m, dst_row, src_row, nb_rows);
}


CML_API void cml_cpy_self_col(MATRIX *m, size_t dst_col, size_t src_col) {
    if (m == NULL || dst_col >= m->cols || src_col >= m->cols) {
        errno = EINVAL;
        return;
    }

    for (size_t i = 0; i < m->rows; ++i) {
        cml_set(m, i, dst_col, cml_get(m, i, src_col));
    }
}


CML_API void cml_cpy_self_col_range(MATRIX *m, size_t dst_col, size_t src_col, size_t nb_cols) {
    if (m == NULL || dst_col + nb_cols > m->cols || src_col + nb_cols >=  m->cols) {
        errno = EINVAL;
        return;
    }
    
    for(size_t i = 0; i < nb_cols; ++i) {
        // void * memcpy ( void * dstination, const void * source, size_t num );
        memcpy(
            (void *)(MATRIX_COL2PTR(m, dst_col + i)),
            (void *)(MATRIX_COL2PTR(m, src_col + i)),
            (m->cols * sizeof(cml_real_t))
        );
    }
}


CML_API void cml_cpy_self_elem(MATRIX *m, size_t dst_row, size_t dst_col, size_t src_row, size_t src_col) {
    if (m == NULL || dst_row >= m->rows || dst_col >= m->cols || src_row >= m->rows || src_col >= m->cols) {
        errno = EINVAL;
        return;
    }

    cml_set(m, dst_row, dst_col, cml_get(m, src_row, src_col));
}


CML_API void cml_cpy_submat(MATRIX *dst, MATRIX *src, size_t dst_row, size_t dst_col, size_t src_row, size_t src_col, size_t nb_rows, size_t nb_cols) {
    /* TODO checks... */

    // for(size_t j = 0; j < nb_cols; ++j) {
    //     // void * memcpy ( void * destination, const void * source, size_t num );
    //     memcpy(
    //         (      void *)(MATRIX_COL2PTR(dst, dst_col + j) + dst_row),
    //         (const void *)(MATRIX_COL2PTR(src, src_col + j) + src_row),
    //         (nb_rows * sizeof(cml_real_t))
    //     );
    // }

    for(size_t j = 0; j < nb_cols; ++j) {
        for(size_t i = 0; i < nb_cols; ++i) {
            cml_set(
                dst,
                dst_row + i,
                dst_col + j,
                cml_get(
                    src,
                    src_row + i,
                    src_col + j
                )
            );
        }
    }
}


CML_API void cml_swap(MATRIX *m1, MATRIX *m2) {
    if (m1 == NULL || m2 == NULL || m1 == m2) {
        errno = EINVAL;
        return;
    }

    MATRIX temp = { m1->rows, m1->cols, m1->data };

    m1->rows = m2->rows;
    m1->cols = m2->cols;
    m1->data = m2->data;

    m2->rows = temp.rows;
    m2->cols = temp.cols;
    m2->data = temp.data;
}


CML_API void cml_swap_row(MATRIX *m1, MATRIX *m2, size_t m1_row, size_t m2_row) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1->cols != m2->cols || m1_row >= m1->rows || m2_row >= m2->rows) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(1, m1->cols);

    if (temp == NULL) {
        return;
    }

    cml_cpy_row(temp, m1, 0, m1_row);
    cml_cpy_row(m1, m2, m1_row, m2_row);
    cml_cpy_row(m2, temp, m2_row, 0);
    cml_free(temp);
}


CML_API void cml_swap_col(MATRIX *m1, MATRIX *m2, size_t m1_col, size_t m2_col) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1->rows != m2->rows || m1_col >= m1->cols || m2_col >= m2->cols) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(m1->rows, 1);

    if (temp == NULL) {
        return;
    }

    cml_cpy_col(temp, m1, 0, m1_col);
    cml_cpy_col(m1, m2, m1_col, m2_col);
    cml_cpy_col(m2, temp, m2_col, 0);
    cml_free(temp);
}


CML_API void cml_swap_elem(MATRIX *m1, MATRIX *m2, size_t m1_row, size_t m1_col, size_t m2_row, size_t m2_col) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1_row >= m1->rows || m1_col >= m1->cols || m2_row >= m2->rows || m2_col >= m2->cols) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(1, 1);

    if (temp == NULL) {
        return;
    }

    cml_cpy_elem(temp, m1, 0, 0, m1_row, m1_col);
    cml_cpy_elem(m1, m2, m1_row, m1_col, m2_row, m2_col);
    cml_cpy_elem(m2, temp, m2_row, m2_col, 0, 0);
    cml_free(temp);
}


CML_API void cml_swap_self_row(MATRIX *m, size_t row_1, size_t row_2) {
    if (m == NULL || row_1 >= m->rows || row_2 >= m->rows) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(1, m->cols);

    if (temp == NULL) {
        return;
    }

    cml_cpy_row(temp, m, 0, row_1);
    cml_cpy_self_row(m, row_1, row_2);
    cml_cpy_row(m, temp, row_2, 0);
    cml_free(temp);
}


CML_API void cml_swap_self_col(MATRIX *m, size_t col_1, size_t col_2) {
    if (m == NULL || col_1 >= m->cols || col_2 >= m->cols) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(m->rows, 1);

    if (temp == NULL) {
        return;
    }

    cml_cpy_col(temp, m, 0, col_1);
    cml_cpy_self_col(m, col_1, col_2);
    cml_cpy_col(m, temp, col_2, 0);
    cml_free(temp);
}


CML_API void cml_swap_self_elem(MATRIX *m, size_t row_1, size_t col_1, size_t row_2, size_t col_2) {
    if (m == NULL || row_1 >= m->rows || col_1 >= m->cols || row_2 >= m->rows || col_2 >= m->cols) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(1, 1);

    if (temp == NULL) {
        return;
    }

    cml_cpy_elem(temp, m, 0, 0, row_1, col_1);
    cml_cpy_self_elem(m, row_1, col_1, row_2, col_2);
    cml_cpy_elem(m, temp, row_2, col_2, 0, 0);
    cml_free(temp);
}


CML_API void cml_ins_row(MATRIX *dst, MATRIX *src, size_t dst_row, size_t src_row) {
    if (dst == NULL || src == NULL || dst == src || dst->cols != src->cols || dst_row > dst->rows || src_row >= src->rows) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dst->data, vec_len(dst->data) + dst->cols);

    if (vec_cap(dst->data) < vec_len(dst->data) + dst->cols) {
        errno = ENOMEM;
        return;
    }

    for (size_t j = dst->cols; j-- > 0; ) {
        vec_insert(dst->data, cml_get(src, src_row, j), dst->rows * j + dst_row);
    }
    dst->rows += 1;
}


CML_API void cml_ins_col(MATRIX *dst, MATRIX *src, size_t dst_col, size_t src_col) {
    if (dst == NULL || src == NULL || dst == src || dst->rows != src->rows || dst_col > dst->cols || src_col >= src->cols) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dst->data, vec_len(dst->data) + dst->rows);
    
    if (vec_cap(dst->data) < vec_len(dst->data) + dst->rows) {
        errno = ENOMEM;
        return;
    }

    for (size_t i = dst->rows; i-- > 0; ) {
        vec_insert(dst->data, cml_get(src, i, src_col), dst->rows * dst_col);
    }
    dst->cols += 1;
}


CML_API void cml_ins_self_row(MATRIX *m, size_t dst_row, size_t src_row) {
    if (m == NULL || dst_row > m->rows || src_row >= m->rows) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(1, m->cols);

    if (temp == NULL) {
        return;
    }

    cml_cpy_row(temp, m, 0, src_row);
    cml_ins_row(m, temp, dst_row, 0);
    cml_free(temp);
}


CML_API void cml_ins_self_col(MATRIX *m, size_t dst_col, size_t src_col) {
    if (m == NULL || dst_col > m->cols || src_col >= m->cols) {
        errno = EINVAL;
        return;
    }

    MATRIX *temp = cml_new(m->rows, 1);

    if (temp == NULL) {
        return;
    }

    cml_cpy_col(temp, m, 0, src_col);
    cml_ins_col(m, temp, dst_col, 0);
    cml_free(temp);
}


CML_API void cml_del_all(MATRIX *m) {
    if (m == NULL) {
        errno = EINVAL;
        return;
    }

    m->rows = 0;
    m->cols = 0;
    vec_clear(m->data);
    vec_shrink(m->data);
}


CML_API void cml_del_row(MATRIX *m, size_t row) {
    if (m == NULL || row >= m->rows) {
        errno = EINVAL;
        return;
    }

    for (size_t j = m->cols; j-- > 0; ) {
        vec_remove(m->data, m->rows * j + row);
    }
    m->rows -= 1;
    vec_shrink(m->data);
}


CML_API void cml_del_col(MATRIX *m, size_t col) {
    if (m == NULL || col >= m->cols) {
        errno = EINVAL;
        return;
    }

    for (size_t i = m->rows; i-- > 0; ) {
        vec_remove(m->data, m->rows * col);
    }
    m->cols -= 1;
    vec_shrink(m->data);
}


CML_API void cml_adjoin_top(MATRIX *dst, MATRIX *src) {
    if (dst == NULL || src == NULL || dst == src || dst->cols != src->cols) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dst->data, vec_len(dst->data) + vec_len(src->data));

    if (vec_cap(dst->data) < vec_len(dst->data) + vec_len(src->data)) {
        errno = ENOMEM;
        return;
    }
    
    for (size_t i = src->rows; i-- > 0; ) {
        cml_ins_row(dst, src, 0, i);
    }
}


CML_API void cml_adjoin_bottom(MATRIX *dst, MATRIX *src) {
    if (dst == NULL || src == NULL || dst == src || dst->cols != src->cols) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dst->data, vec_len(dst->data) + vec_len(src->data));

    if (vec_cap(dst->data) < vec_len(dst->data) + vec_len(src->data)) {
        errno = ENOMEM;
        return;
    }

    for (size_t i = 0; i < src->rows; ++i) {
        cml_ins_row(dst, src, dst->rows, i);
    }
}


CML_API void cml_adjoin_left(MATRIX *dst, MATRIX *src) {
    if (dst == NULL || src == NULL || dst == src || dst->rows != src->rows) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dst->data, vec_len(dst->data) + vec_len(src->data));

    if (vec_cap(dst->data) < vec_len(dst->data) + vec_len(src->data)) {
        errno = ENOMEM;
        return;
    }

    for (size_t j = src->cols; j-- > 0; ) {
        cml_ins_col(dst, src, 0, j);
    }
}


CML_API void cml_adjoin_right(MATRIX *dst, MATRIX *src) {
    if (dst == NULL || src == NULL || dst == src || dst->rows != src->rows) {
        errno = EINVAL;
        return;
    }

    vec_reserve(dst->data, vec_len(dst->data) + vec_len(src->data));

    if (vec_cap(dst->data) < vec_len(dst->data) + vec_len(src->data)) {
        errno = ENOMEM;
        return;
    }
    
    for (size_t j = 0; j < src->cols; ++j) {
        cml_ins_col(dst, src, dst->cols, j);
    }
}


CML_API void cml_add_const(MATRIX *m, cml_real_t data, MATRIX *opt) {
    if (m == NULL || m == opt || m->rows == 0) {
        errno = EINVAL;
        return;
    }

    if (opt != NULL) {
        cml_cpy(opt, m);
        m = opt;
    } 
    for (size_t i = 0; i < vec_len(m->data); ++i) {
        m->data[i] += data;
    }
}


CML_API void cml_mul_const(MATRIX *m, cml_real_t data, MATRIX *opt) {
    if (m == NULL || m == opt || m->rows == 0) {
        errno = EINVAL;
        return;
    }

    if (opt != NULL) {
        cml_cpy(opt, m);
        m = opt;
    } 
    for (size_t i = 0; i < vec_len(m->data); ++i) {
        m->data[i] *= data;
    }
}


CML_API void cml_add(MATRIX *m1, MATRIX *m2, MATRIX *opt) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1 == opt || m2 == opt || m1->rows != m2->rows || m1->cols != m2->cols || m1->rows == 0) {
        errno = EINVAL;
        return;
    }

    if (opt != NULL) {
        cml_cpy(opt, m1);
        m1 = opt;
    }
    for (size_t i = 0; i < vec_len(m1->data); ++i) {
        m1->data[i] += m2->data[i];
    }
}


CML_API void cml_sub(MATRIX *m1, MATRIX *m2, MATRIX *opt) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1 == opt || m2 == opt || m1->rows != m2->rows || m1->cols != m2->cols || m1->rows == 0) {
        errno = EINVAL;
        return;
    }
    
    if (opt != NULL) {
        cml_cpy(opt, m1);
        m1 = opt;
    }
    for (size_t i = 0; i < vec_len(m1->data); ++i) {
        m1->data[i] -= m2->data[i];
    }
}


CML_API void cml_mul(MATRIX *m1, MATRIX *m2, MATRIX *opt) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1 == opt || m2 == opt || m1->cols != m2->rows || m1->rows == 0) {
        errno = EINVAL;
        return;
    }

    MATRIX *result = cml_new(m1->rows, m2->cols);

    if (result == NULL) {
        return;
    } 

    #ifndef CML_NO_DEPENDENCIES
    /* TODO - handle single vs. double precision */
    char notrans[2] = "N";
    int m = (int) m1->rows;
    int n = (int) m2->cols;
    int k = (int) m1->cols;
    cml_real_t alpha = 1;
    cml_real_t beta = 0;

    dgemm_(notrans, notrans, &m, &n, &k, &alpha, m1->data, &m, m2->data, &k, &beta, result->data, &m);
    #else
    for (size_t i = 0; i < m1->rows; ++i) {
        for (size_t j = 0; j < m2->cols; ++j) {
            cml_real_t sum = 0;
            for (size_t k = 0; k < m1->cols; ++k) {
                sum += cml_get(m1, i, k) * cml_get(m2, k, j);
            }
            cml_set(result, i, j, sum);
        }
    }
    #endif
    
    if (opt != NULL) {
        cml_cpy(opt, result);
    } else {
        cml_cpy(m1, result);
    }
    cml_free(result);
}


CML_API void cml_mul_elem(MATRIX *m1, MATRIX *m2, MATRIX *opt) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1 == opt || m2 == opt || m1->rows != m2->rows || m1->cols != m2->cols || m1->rows == 0) {
        errno = EINVAL;
        return;
    }
    
    if (opt != NULL) {
        cml_cpy(opt, m1);
        m1 = opt;
    }
    for (size_t i = 0; i < vec_len(m1->data); ++i) {
        m1->data[i] *= m2->data[i];
    }
}


CML_API void cml_div_elem(MATRIX *m1, MATRIX *m2, MATRIX *opt) {
    if (m1 == NULL || m2 == NULL || m1 == m2 || m1 == opt || m2 == opt || m1->rows != m2->rows || m1->cols != m2->cols || m1->rows == 0) {
        errno = EINVAL;
        return;
    }
    
    if (opt != NULL) {
        cml_cpy(opt, m1);
        m1 = opt;
    }
    for (size_t i = 0; i < vec_len(m1->data); ++i) {
        m1->data[i] /= m2->data[i];
    }
}

/* add two channels in place */
CML_API void cml_add_chan_m1_m2(MATRIX *m1, size_t m1ch, MATRIX *m2, size_t m2ch) {
    if (m1 == NULL || m2 == NULL || m1->rows != m2->rows || m1->rows == 0 || m1ch >= m1->rows || m2ch >= m2->rows) {
        errno = EINVAL;
        return;
    }

    /* opt can be equal to m1 or m2, but must be specified by caller */
    cml_real_t *m1ptr = MATRIX_CHAN2PTR(m1, m1ch);
    cml_real_t *m2ptr = MATRIX_CHAN2PTR(m2, m2ch);
    for (size_t i = 0; i < m1->rows; ++i) {
        *(m1ptr++) += *(m2ptr++);
    }
}


CML_API cml_real_t cml_sum(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return 0;
    }

    cml_real_t sum = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            sum += cml_get(m, i, j);
        }
    }
    return sum;
}


CML_API void cml_sum_dim(MATRIX *m, size_t dim, MATRIX *opt) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return;
    }
    /* TODO assert (or errno) size checks */

    if(0 == dim) {
        for (size_t j = 0; j < m->cols; ++j) {
            cml_real_t sum = 0;
            for (size_t i = 0; i < m->rows; ++i) {
                sum += cml_get(m, i, j);
            }
        }

    }
    else if(1 == dim) {
        for (size_t i = 0; i < m->rows; ++i) {
            cml_real_t sum = 0;
            for (size_t j = 0; j < m->cols; ++j) {
                sum += cml_get(m, i, j);
            }
            cml_set(opt, i, 0, sum);
        }
    }
}


CML_API cml_real_t cml_mean(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return 0;
    }
    return cml_sum(m) / (m->rows * m->cols);
}


CML_API void cml_mean_dim(MATRIX *m, size_t dim, MATRIX *opt) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return;
    }
    /* TODO assert (or errno) size checks */

    cml_sum_dim(m, dim, opt);
    if(0 == dim) {
        cml_mul_const(opt, (cml_real_t)1/(cml_real_t)m->rows, NULL);
    }
    else if(1 == dim) {
        cml_mul_const(opt, (cml_real_t)1/(cml_real_t)m->cols, NULL);
    }
}


CML_API cml_real_t cml_rms(MATRIX *m){
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return 0;
    }

    cml_real_accum_t tmp = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            tmp += (cml_real_accum_t)cml_real_square(cml_get(m, i, j));
        }
    }
    return cml_real_sqrt((cml_real_t)(tmp / (cml_real_accum_t)(m->rows * m->cols)));
}


// CML_API  void          cml_rms_dim(MATRIX *m, MATRIX *opt, size_t dim);


CML_API void cml_abs(MATRIX *m, MATRIX *opt) {
    cml_elm_real_func0(m, cml_real_abs, opt);
}


CML_API void cml_pow_mat2const(MATRIX *m, cml_real_t r1, MATRIX *opt) {
    cml_elm_real_func1(m, cml_real_pow, r1, opt);
}


CML_API void cml_pow_const2mat(cml_real_t r1, MATRIX *m, MATRIX *opt) {
    cml_elm_real_1func(r1, cml_real_pow, m, opt);
}


CML_API void cml_log10(MATRIX *m, MATRIX *opt) {
    cml_elm_real_func0(m, cml_real_log10, opt);
}


CML_API void cml_10_log10_abs(MATRIX *m, MATRIX *opt) {
    cml_elm_real_func0(m, cml_real_10log10abs, opt);
}


CML_API void cml_20_log10_abs(MATRIX *m, MATRIX *opt) {
    cml_elm_real_func0(m, cml_real_20log10abs, opt);
}


CML_API void cml_clip(MATRIX *m, cml_real_t lolim, cml_real_t hilim, MATRIX *opt) {
    cml_elm_real_func2(m, cml_real_clip, lolim, hilim, opt);
}


CML_API void cml_gt(MATRIX *m, cml_real_t r1, MATRIX *opt) {
    cml_elm_real_func1(m, cml_real_gt, r1, opt);
}


CML_API void cml_gte(MATRIX *m, cml_real_t r1, MATRIX *opt) {
    cml_elm_real_func1(m, cml_real_gte, r1, opt);
}


CML_API void cml_lt(MATRIX *m, cml_real_t r1, MATRIX *opt) {
    cml_elm_real_func1(m, cml_real_lt, r1, opt);
}


CML_API void cml_lte(MATRIX *m, cml_real_t r1, MATRIX *opt) {
    cml_elm_real_func1(m, cml_real_lte, r1, opt);
}


CML_API cml_real_t cml_min(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return 0;
    }

    cml_real_t min = cml_get(m, 0, 0);
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            if (min > cml_get(m, i, j)) {
                min = cml_get(m, i, j);
            }
        }
    }

    return min;
}


CML_API void cml_min_dim(MATRIX *m, size_t dim, MATRIX *opt) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return;
    }

    if(dim==0) { /* min per 0th dimension */
        /* TODO assert (or errno) size checks */
        for (size_t j = 0; j < m->cols; ++j) {
            cml_real_t min = cml_get(m, 0, j);
            for (size_t i = 0; i < m->cols; ++i) {
                if (min > cml_get(m, i, j)) {
                    min = cml_get(m, i, j);
                }
            }
            cml_set(opt, j, 0, min);
        }
    }
    else if(dim==1) { /* min per 0th dimension */
        /* TODO assert (or errno) size checks */
        for (size_t i = 0; i < m->rows; ++i) {
            cml_real_t min = cml_get(m, i, 0);
            for (size_t j = 0; j < m->cols; ++j) {
                if (min > cml_get(m, i, j)) {
                    min = cml_get(m, i, j);
                }
            }
            cml_set(opt, i, 0, min);
        }
    }
}



CML_API cml_real_t cml_max(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return 0;
    }

    cml_real_t max = cml_get(m, 0, 0);
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            if (max < cml_get(m, i, j)) {
                max = cml_get(m, i, j);
            }
        }
    }

    return max;
}


CML_API void cml_max_dim(MATRIX *m, size_t dim, MATRIX *opt) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return;
    }

    if(dim==0) { /* max per 0th dimension */
        /* TODO assert (or errno) size checks */
        for (size_t j = 0; j < m->cols; ++j) {
            cml_real_t max = cml_get(m, 0, j);
            for (size_t i = 0; i < m->cols; ++i) {
                if (max < cml_get(m, i, j)) {
                    max = cml_get(m, i, j);
                }
            }
            cml_set(opt, j, 0, max);
        }
    }
    else if(dim==1) { /* max per 0th dimension */
        /* TODO assert (or errno) size checks */
        for (size_t i = 0; i < m->rows; ++i) {
            cml_real_t max = cml_get(m, i, 0);
            for (size_t j = 0; j < m->cols; ++j) {
                if (max < cml_get(m, i, j)) {
                    max = cml_get(m, i, j);
                }
            }
            cml_set(opt, i, 0, max);
        }
    }
    return;
}


CML_API bool cml_is_zero(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return false;
    }

    for (size_t i = 0; i < vec_len(m->data); ++i) {
        if (m->data[i] != 0) {
            return false;
        }
    }

    return true;
}


// CML_API bool cml_is_zero_dim(MATRIX *m, size_t dim, MATRIX *opt);


CML_API bool cml_is_pos(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return false;
    }

    for (size_t i = 0; i < vec_len(m->data); ++i) {
        if (m->data[i] <= 0) {
            return false;
        }
    }

    return true;
}


// CML_API bool cml_is_pos_dim(MATRIX *m, size_t dim, MATRIX *opt) {}


CML_API bool cml_is_nonneg(MATRIX *m) {
    if (m == NULL || m->rows == 0) {
        errno = EINVAL;
        return false;
    }

    for (size_t i = 0; i < vec_len(m->data); ++i) {
        if (m->data[i] < 0) {
            return false;
        }
    }

    return true;
}


// CML_API bool cml_is_pos_nonneg_dim(MATRIX *m, size_t dim, MATRIX *opt) {}


CML_API bool cml_is_equal(MATRIX *m1, MATRIX *m2) {
    if (m1 == NULL || m2 == NULL || m1 == m2) {
        errno = EINVAL;
        return false;
    }

    if (m1->rows != m2->rows || m1->cols != m2->cols) {
        return false;
    }
    for (size_t i = 0; i < vec_len(m1->data); ++i) {
        if (m1->data[i] != m2->data[i]) {
            return false;
        }
    }
    
    return true;
}


CML_API void cml_transpose(MATRIX *m, MATRIX *opt) {
    if (m == NULL || m == opt) {
        errno = EINVAL;
        return;
    }

    if (opt != NULL) {
        cml_cpy(opt, m);
        m = opt;
    }
    MATRIX *trans = cml_new(m->cols, m->rows);
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            cml_set(trans, j, i, cml_get(m, i, j));
        }
    }
    cml_cpy(m, trans);
    cml_free(trans);
}


CML_API void cml_normalize(MATRIX *m, MATRIX *opt) {
    if (m == NULL || m == opt || m->rows == 0) {
        errno = EINVAL;
        return;
    }

    if (opt != NULL) {
        cml_cpy(opt, m);
        m = opt;
    }
    cml_real_t min = cml_min(m);
    cml_real_t max = cml_max(m);
    for (size_t i = 0; i < vec_len(m->data); ++i) {
        m->data[i] = (m->data[i] - min) / (max - min);
    }
}


CML_API void cml_normalize_dim(MATRIX *m, size_t dim, MATRIX *opt) {
    if (m == NULL || m == opt || m->rows == 0 || (dim!=0 && dim!=1)) {
        errno = EINVAL;
        return;
    }

    size_t mmlen = dim ? m->rows : m->cols;
    /* TODO assert (or errno) dim sizes */
    MATRIX *vecmin = cml_new(mmlen, 1);
    MATRIX *vecmax = cml_new(mmlen, 1);

    if (opt == NULL) {
        opt = m;
    }

    /* store min and max in temp vectors */
    cml_min_dim(m, dim, vecmin);
    cml_max_dim(m, dim, vecmax);

    if(0 == dim) {
        for(size_t j=0; j<m->cols; j++) {
            cml_real_t min=cml_get(vecmin, j, 0), max=cml_get(vecmax, j, 0);
            for(size_t i=0; i<m->rows; i++) {
                cml_set(opt, i, j, (cml_get(m, i, j) - min) / (max - min));
            }
        }
    }
    else if(1 == dim) {
        for(size_t i=0; i<m->rows; i++) {
            cml_real_t min=cml_get(vecmin, i, 0), max=cml_get(vecmax, i, 0);
            for(size_t j=0; j<m->cols; j++) {
                cml_set(opt, i, j, (cml_get(m, i, j) - min) / (max - min));
            }
        }
    }

    /* free temporary vectors */
    cml_free(vecmin);
    cml_free(vecmax);
}

/*    SIGNAL PROCESSING (ETC)    */
CML_API void __biquad_proc(cml_real_t *bq, cml_real_t *in, cml_real_t *out, size_t len) {
    cml_real_t d0;
    cml_real_t *w = bq + 5;
    for (size_t fix=0; fix<len; fix++) {
        d0 =   in[fix]      \
             - bq[3] * w[0] \
             - bq[4] * w[1] ;
        out[fix] =  bq[0] *   d0 \
                   + bq[1] * w[0] \
                   + bq[2] * w[1] ;
        w[1] = w[0];
        w[0] = d0;
    }
}

CML_API void __sos_proc(cml_real_t *sos, size_t nb_bq, cml_real_t *in, size_t len, cml_real_t *out, cml_real_t *scratch) {
    bool free_scratch = false;
    size_t ppix = 0;
    cml_real_t *buf[2], *src, *dst;

    /* determine if we should malloc a scratch buffer */
    if(NULL == scratch) {
        scratch = (cml_real_t *)calloc(len, sizeof(cml_real_t));
        free_scratch = true;
    }

    /* set up ping-pong buffers between output and scratch */
    buf[0] = out;
    buf[1] = scratch;

    /* run biquads[0] from input to output */
    src=in; dst=out;

    for (size_t bqix = 0; bqix < nb_bq; ++bqix) {
        __biquad_proc(sos + (7*bqix), src, dst, len);
        src=buf[ppix&1]; dst=buf[(ppix+1)&1];
        ppix++;
    }

    /* if necessary, copy the final result in scratch to output */
    if(dst != out && src == out ) {
        // void * memcpy ( void * dstination, const void * source, size_t num );
        memcpy((void *)out, (void *)scratch, len * sizeof(cml_real_t));
    }

    /* free scratch if it was allocated in this function */
    if( free_scratch ) {
        free(scratch);
    }
}

/* TODO - just split up these functions for better readability */
CML_API void cml0_sos_proc(MATRIX *sos, MATRIX *in, MATRIX *scratch, MATRIX *out, size_t in_ch, size_t out_ch) {
    bool free_scratch = false;
    if (sos == NULL || in == NULL || in->rows == 0) {
        errno = EINVAL;
        return;
    }
    if(out == NULL) {
        out = in;
    }

    if(DIFF_DIM_SZ(in, 0, out, 0)) {
        errno = EINVAL;
        return;
    }

    if(NULL == scratch) {
        scratch = cml_new(in->rows, 1);
        free_scratch = true;
    }

    size_t nb_bq = sos->rows / 7;

    /* 1 to 1 */
    if (in_ch < in->cols && out_ch < out->cols) {
        __sos_proc(
            MATRIX_CHAN2PTR(sos, 0), nb_bq,
            MATRIX_CHAN2PTR(in, in_ch), in->rows,
            scratch->data,
            MATRIX_CHAN2PTR(out, out_ch)
        );
    }

    /* ALL to ALL */
    else if (in_ch == MATRIX_ALL_CH && out_ch == MATRIX_ALL_CH) {
        if(DIFF_DIM_SZ(in, 1, out, 1)) {
            errno = EINVAL;
            return;
        }
        for(size_t j=0; j<in->cols; j++) {
            __sos_proc(
                MATRIX_CHAN2PTR(sos, j), nb_bq,
                MATRIX_CHAN2PTR(in, j), in->rows,
                scratch->data,
                MATRIX_CHAN2PTR(out, j)
            );
        }
    }

    /* ALL to 1 */
    else if (in_ch == MATRIX_ALL_CH && out_ch < out->cols) {
        /* in == out is invalid for this case of 'ALL to 1' */
        if (in == out ) {
            errno = EINVAL;
            return;
        }
        MATRIX *tmpmat = cml_new(in->rows, 1);
        cml_set_col(out, out_ch, 0); /* zero final dstination */
        for (size_t j = 0; j < in->cols; ++j) {
            __sos_proc(
                MATRIX_CHAN2PTR(sos, j), nb_bq,
                MATRIX_CHAN2PTR(in, j), in->rows,
                scratch->data,
                MATRIX_CHAN2PTR(tmpmat, 0)
            );
            /* sum result from tmpmat in to out_ch of out */
            cml_add_chan_m1_m2(out, out_ch, tmpmat, 0);
        }
        cml_free(tmpmat);
    }

    /* 1 to ALL */
    else if (in_ch == MATRIX_ALL_CH && out_ch == MATRIX_ALL_CH) {
        /* filter in to channel 0 of out */
        __sos_proc(
            MATRIX_CHAN2PTR(sos, 0), nb_bq,
            MATRIX_CHAN2PTR(in, in_ch), in->rows,
            scratch->data,
            MATRIX_CHAN2PTR(out, 0)
        );

        /* copy channel 0 of out to all other channels of out */
        for (size_t j = 1; j < out->cols; ++j) {
            cml_cpy_self_col(out, j, 0);
        }
    }

    if(free_scratch) {
        cml_free(scratch);
    }
}


CML_API cml_real_t __dot_prod(cml_real_t *p1, cml_real_t *p2, size_t len, bool flip) {
    cml_real_t out = 0;
    int incr = 1;
    if (flip) {
        incr = -1;
        p2 += (len - 1);
    }
    for (size_t i = 0; i < len; ++i) {
        out += (*(p1++)) * (*p2);
        p2 += incr;
    }
    return out;
}


CML_API void cml0_fir_proc(MATRIX *coef, MATRIX *state, MATRIX *in, MATRIX *out, size_t len) {
    if (NULL == coef || NULL == state || NULL == in || NULL == out || in->cols != out->cols || \
            0 == coef->rows || coef->rows + len - 1 > state->rows || len > out->rows) {
        errno = EINVAL;
        return;
    }

    /* copy input to state matrix where previous input is at the beginning */
    cml_cpy_row_range(state, in, coef->rows - 1, 0, len);

    for(size_t j = 0; j < in->cols; ++j) {
        /* filter each channel */
        cml_real_t *cptr = (j < coef->cols) ? MATRIX_CHAN2PTR(coef, j) : MATRIX_CHAN2PTR(coef, 0);
        for(size_t i = 0; i < len; ++i) {
            cml_set(out,
                i, j,
                __dot_prod(
                    MATRIX_CHAN2PTR(state, j) + i, // data pointer
                    cptr,                          // coef pointer
                    coef->rows,                    // len of dot_prod / conv
                    true                           // flip to do conv instead of dot_prod
                )
            );
        }
    }

    /* copy the last part of state to the beginning for the next call */
    cml_cpy_row_range(state, state, 0, len - coef->rows + 1, len - 1);
}


#ifndef CML_NO_DEPENDENCIES


/* Returns false if matrix m doesn't have inverse or error occurred. Otherwise returns true.
 */
CML_API bool cml_inverse(MATRIX *m, MATRIX *opt) {
    if (m == NULL || m == opt || m->rows != m->cols || m->rows == 0) {
        errno = EINVAL;
        return false;
    }

    int old_errno = errno;
    errno = 0;
    MATRIX *backup = cml_new(0, 0);
    cml_cpy(backup, m);
    if (opt != NULL) {
        cml_cpy(opt, m);
        m = opt;
    } 

    if (errno) {
        cml_free(backup);
        return false;
    } else {
        errno = old_errno;
    }

    int n = (int) m->rows;
    int ipiv[n];
    int lwork = n * n;
    cml_real_t work[lwork];
    int info;
    dgetrf_(&n, &n, m->data, &n, ipiv, &info);
    dgetri_(&n, m->data, &n, ipiv, work, &lwork, &info);

    if (info < 0) {
        cml_cpy(m, backup);
        cml_free(backup);
        errno = EDOM;
        return false;
    } else if (info > 0) {
        cml_cpy(m, backup);
        cml_free(backup);
        return false;
    }

    cml_free(backup);

    return true;
}


/* Solves the system of equations 'AX = B' for X and stores it in B (NOT A, since B already has desired dimensions), or
 * 'opt' if 'opt != NULL'. Returns false if unsolvable or error occurred. Otherwise returns true.
 */
CML_API bool cml_sys_equ(MATRIX *A, MATRIX *B, MATRIX *opt) {
    if (A == NULL || B == NULL || A == B || A == opt || B == opt || A->rows != A->cols || A->cols != B->rows || A->rows == 0) {
        errno = EINVAL;
        return false;
    }

    int old_errno = errno;
    errno = 0;
    MATRIX *A_copy = cml_new(0, 0);
    cml_cpy(A_copy, A);
    A = A_copy;
    MATRIX *backup = cml_new(0, 0);
    cml_cpy(backup, B);
    if (opt != NULL) {
        cml_cpy(opt, B);
        B = opt;
    }

    if (errno) {
        cml_free(A);
        cml_free(backup);
        return false;
    } else {
        errno = old_errno;
    }

    int n = (int) A->rows;
    int nrhs = (int) B->cols;
    int ldb = (int) B->rows;
    int ipiv[n];
    int info;
    dgesv_(&n, &nrhs, A->data, &n, ipiv, B->data, &ldb, &info);

    if (info < 0) {
        cml_cpy(B, backup);
        cml_free(A);
        cml_free(backup);
        errno = EDOM;
        return false;
    } else if (info > 0) {
        cml_cpy(B, backup);
        cml_free(A);
        cml_free(backup);
        return false;
    }

    cml_free(A);
    cml_free(backup);

    return true;
}


#endif


#endif


#endif 
