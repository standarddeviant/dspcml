
#ifndef CML_SAMPLE_TYPE_H
#define CML_SAMPLE_TYPE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

/*    SAMPLE TYPE    */
#ifndef cml_real_t
#define cml_real_t float
#define cml_sample_pow powf
#define cml_sample_log10 log10f
#define cml_sample_square(x) (x*x)
#define cml_sample_sqrt sqrtf
#define cml_sample_abs fabsf
#endif

#ifdef __cplusplus
}
#endif

#endif // CML_SAMPLE_TYPE_H
