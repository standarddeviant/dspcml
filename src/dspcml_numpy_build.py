# file "dspcml_numpy_build.py"

# MIT License
# Copyright (c) [2019] [David Crist]

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef(
"""
typedef float cml_real_t;
size_t maybe_print_errno(char *cptr);
size_t print_matrix(unsigned long ulptr, char *cptr);
unsigned long zeros(size_t rows, size_t cols);
unsigned long arange(cml_real_t start, cml_real_t stop, cml_real_t step);
unsigned long tile(unsigned long ulptr, size_t row_copies, size_t col_copies);
unsigned long new_sos(size_t nb_bq, size_t nb_ch, cml_real_t *coeffs);
void cpy_submat(unsigned long ulptr_dst, unsigned long ulptr_src, size_t dst_row, size_t dst_col, size_t src_row, size_t src_col, size_t nb_rows, size_t nb_cols);
size_t rows(unsigned long ulptr);
size_t cols(unsigned long ulptr);
size_t numel(unsigned long ulptr);
void to_numpy(unsigned long ulptr, cml_real_t *np_data);
unsigned long to_dspcml(cml_real_t *in_data, size_t rows, size_t cols);
void __sin__(unsigned long ulptr, unsigned long ulptr_opt);
void __cos__(unsigned long ulptr, unsigned long ulptr_opt);
void reshape(unsigned long ulptr, size_t new_rows, size_t new_cols);
"""
)

ffibuilder.set_source("_dspcml_numpy",  # name of the output C extension
    """
    #include "dspcml_element_type.h"
    #include "dspcml_numpy.h"
    //
    """,
    sources=["dspcml_numpy.c"],   # includes pi.c as additional sources
    libraries=["m"])    # on Unix, link with the math library

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)