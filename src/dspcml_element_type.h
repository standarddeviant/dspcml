
#ifndef CML_ELEMENT_TYPE_H
#define CML_ELEMENT_TYPE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

/*    SAMPLE TYPES    */
#ifndef cml_real_t
#define cml_real_t float
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

cml_real_t cml_real_10log10abs(cml_real_t x);
cml_real_t cml_real_20log10abs(cml_real_t x);
cml_real_t cml_real_clip(cml_real_t x, cml_real_t lolim, cml_real_t hilim);


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
#endif

#ifdef __cplusplus
}
#endif

#endif // CML_ELEMENT_TYPE_H
