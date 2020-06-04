# NOTE: This repository should not be used for anything serious unless and until there are mature tests for correctness.

# DSPCML
DSPCML is an acronym for Digital Signal Processing C Matrix Library. It is an adaption of the C Matrix Library (CML) referenced here: https://github.com/MichaelJWelsh/cml

DSPCML aims to add common signal processing functions to the original CML implementation. The primary goal is to achieve portability with ability to support optimized signal processing functions for different compilers and toolchains via preprocessor options, i.e. `#ifdef XYZ`. One example of this is when using the ESP-IDF for the ESP32 chip, the preprocessor check `#if ESP_PLATFORM` can be used to determine if building for the ESP32 or ESP32-S2 chips. In this cast, the optimized functions to be used are contained in `esp-dsp` : https://github.com/espressif/esp-dsp

## Notable changes from CML to DSPCML
* The elementary data type of CML was `double`. In DSPCML this is changed to `cml_real_t` with the default as `float`. This type alias, along with per-sample math functions, is defined `cml_real_type.h`
* Support for a complex data type, `cml_cpx_t` is being added. This is simply simply two instances of `cml_real_t` that represents the real and imaginary parts of a complex number.
* The original CML implementation used the `cml` prefix in function names. DSPCML adds a few more for the following reasons
    - `cml0` - these functions are explicitly designed to do work along dimension 0, typically as 'channels' of 'samples' in the context of signal processing. This fuzzy distinction is where DSP meets CML. Some of these functions exploit the fact that dimension 0 of a CML matrix is laid out in contiguous memory.
    - `cmlz` - these functions are designed to work on complex matrices, i.e. instances of `ZMATRIX`
    - `cmlz0` - these functions are designed to work on the 0th dimension of complex matrices, i.e. instances of `ZMATRIX`

## Signal Processing Function Wishlist
DSPCML aims to be a strict superset of the functions in CML.
The signal processing function wish-list is maintained below:
#### Per-Element (output same dims as input)
- [x] 10 * LOG10(ABS)
- [x] 20 * LOG10(ABS)
- [x] ABS
- [x] CLIP
- [x] `>`  (GREATER THAN)
- [x] `>=` (GREATER THAN OR EQUAL)
- [x] `<`  (LESS THAN)
- [x] `<=` (LESS THAN OR EQUAL)
- [x] SIN
- [x] COS

#### Per-Dimension (output reduced dims from input)
- [x] MAX
- [x] MEAN
- [x] MIN
- [x] RMS
- [x] SUM

#### Per-Dimension-0 (output same columns as input)
- [x] GENERATE LIN
- [x] GENERATE SIN
- [x] GENERATE COS
- [x] SOS FILTER (cascaded biquads)
- [ ] XMEAN (exclusive mean)
- [ ] FIR FILTER
- [ ] FFT (Real input)
- [ ] FFT (Complex input)
- [ ] DTFT (Real input)
- [ ] DTFT (Complex input)


## Optimized Function Wishtable
DSPCML aims to have optimized versions of functions for different hardware processors.
The optimized function wish-table is maintained below:
Hardware              | Biquad / SOS       | FFT / IFFT     | TBD               |
:-------------------- | :----------------- | :------------- | :---------------- |
ARM-Cortex-A          |                    |                |                   |
ARM-Cortex-M          |                    |                |                   |
ESP32 (Xtensa LX6)    |                    |                |                   |
ESP32-S2 (Xtensa LX7) |                    |                |                   |
X-86                  |                    |                |                   |
X-86-64               |                    |                |                   |


The original CML README file is copied below as DSPCML aims to be a superset of CML.

# CML
CML is a (fully C++ compatible) header-only C matrix library designed with a focus on portability, simplicity, and efficiency written under the POSIX standard. Several computationally-intense functions require/are enhanced by BLAS/LAPACK. The user is able to define optional flags prior to including this library in source code to do things such as removing BLAS/LAPACK as a dependency (while simultaneously removing several functions from API), or giving all functions static storage. CML operates on type 'MATRIX' for all processes. 'MATRIX' is column major so the user can easily interface with the most popular linear-algebra-related API's. 'errno' is used for error handling extensively. 'errno' is set only in the presence of an error. The user is free to disregard 'errno' if an error occurred, but the state of the resultant matrix is undefined. 


## Usage
CML can be used in either header-only mode or implementation mode. The header-only mode is used by default when included and allows including this header in other headers and does not contain the actual implementation. 

The implementation mode requires the user to define the preprocessor macro 'CML_IMPLEMENTATION' in one .c/.cpp file before ```#include```ing this file, e.g.:
 ```C		
#define CML_IMPLEMENTATION
#include "cml.h"
```
IMPORTANT: Everytime you include "cml.h", you have to define the same optional flags. Not doing so will lead to compile-time errors.


## Example
```C
#include <stdlib.h>
#include <stdio.h>

#define CML_IMPLEMENTATION
#include "cml.h"

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

int main(void) {
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
```
