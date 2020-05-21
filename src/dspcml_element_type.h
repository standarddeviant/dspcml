
#ifndef CML_ELEMENT_TYPE_H
#define CML_ELEMENT_TYPE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

/*    REAL SAMPLE TYPES    */


#ifndef cml_real_t
typedef float cml_real_t;
// func pointer types for different numbers of extra args
typedef cml_real_t (*cml_real_func0_t)(cml_real_t);
typedef cml_real_t (*cml_real_func1_t)(cml_real_t, cml_real_t);
typedef cml_real_t (*cml_real_func2_t)(cml_real_t, cml_real_t, cml_real_t);
#define cml_real_pow powf
#define cml_real_log10 log10f
#define cml_real_square(x) (x*x)
#define cml_real_sqrt sqrtf
#define cml_real_abs fabsf
#define cml_real_sin sinf
#define cml_real_cos sinf
#endif


#ifndef cml_real_accum_t
/* can be same or diff based on float vs. fixed, etc. */
#define cml_real_accum_t float
#endif


/*    COMPLEX SAMPLE TYPES    */


#ifndef cml_cpx_t
typedef struct cml_cpx_t {
    cml_real_t re;
    cml_real_t im;
} cml_cpx_t;
#endif


#ifndef cml_cpx_accum_t
/* can be same or diff based on float vs. fixed, etc. */
#define cml_cpx_accum_t cml_cpx_t
#endif


/*    REAL SAMPLE FUNCTION    */
cml_real_t cml_real_10log10abs(cml_real_t x);
cml_real_t cml_real_20log10abs(cml_real_t x);
cml_real_t cml_real_clip(cml_real_t x, cml_real_t lolim, cml_real_t hilim);
cml_real_t cml_real_gt(cml_real_t x, cml_real_t thresh);
cml_real_t cml_real_gte(cml_real_t x, cml_real_t thresh);
cml_real_t cml_real_lt(cml_real_t x, cml_real_t thresh);
cml_real_t cml_real_lte(cml_real_t x, cml_real_t thresh);


/* COMPLEX SAMPLE FUNCTION */
cml_real_t cml_cpx_abs(cml_cpx_t x);


#ifdef CML_IMPLEMENTATION
cml_real_t cml_real_10log10abs(cml_real_t x){
    return 10*cml_real_log10(cml_real_abs(x));
}
cml_real_t cml_real_20log10abs(cml_real_t x){
    return 20*cml_real_log10(cml_real_abs(x));
}
cml_real_t cml_real_clip(cml_real_t x, cml_real_t lolim, cml_real_t hilim) {
    return (x < lolim) ? lolim : ((x > hilim) ? hilim : x);
}
cml_real_t cml_real_gt(cml_real_t x, cml_real_t thresh){
    return (x > thresh) ? 1 : 0;
}
cml_real_t cml_real_gte(cml_real_t x, cml_real_t thresh){
    return (x >= thresh) ? 1 : 0;
}
cml_real_t cml_real_lt(cml_real_t x, cml_real_t thresh){
    return (x < thresh) ? 1 : 0;
}
cml_real_t cml_real_lte(cml_real_t x, cml_real_t thresh){
    return (x <= thresh) ? 1 : 0;
}
#endif

#ifdef __cplusplus
}
#endif

#endif // CML_ELEMENT_TYPE_H
