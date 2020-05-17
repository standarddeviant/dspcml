
#ifndef CML_ELEMENT_TYPE_H
#define CML_ELEMENT_TYPE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

/*    SAMPLE TYPES    */
#ifndef cml_real_t
#define cml_real_t float
#define cml_real_pow powf
#define cml_real_log10 log10f
#define cml_real_square(x) (x*x)
#define cml_real_sqrt sqrtf
#define cml_real_abs fabsf
#endif

#ifdef __cplusplus
}
#endif

#endif // CML_ELEMENT_TYPE_H
