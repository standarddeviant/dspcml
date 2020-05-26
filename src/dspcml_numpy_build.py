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
int print_matrix(unsigned long ulptr, char *cptr);
unsigned long zeros(size_t rows, size_t cols);
unsigned long rows(unsigned long ulptr);
unsigned long cols(unsigned long ulptr);
unsigned long numel(unsigned long ulptr);
void to_numpy(unsigned long ulptr, float *np_data);
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