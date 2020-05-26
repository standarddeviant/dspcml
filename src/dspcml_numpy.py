import numpy as np
from _dspcml_numpy import ffi, lib

py_real_str = "float32"
c_real_str = "float"


def zeros(rows, cols):
    # return DSPCML MATRIX pointer as an unsigned long for simplicity
    return lib.zeros(rows, cols)


def __numel__(ulptr):
    # return number of elements in DSPCML MATRIX
    return lib.numel(ulptr)


def __rows__(ulptr):
    return lib.rows(ulptr)


def __cols__(ulptr):
    return lib.cols(ulptr)


def print_matrix(ulptr):
    strarg = ffi.new("char[4096]")
    lib.print_matrix(ulptr, strarg)
    print(ffi.string(strarg).decode())


def to_numpy(ulptr):
    # Use numpy/ctypes trick to make valid C pointer from numpy array
    # https://stackoverflow.com/questions/16276268/how-to-pass-a-numpy-array-into-a-cffi-function-and-how-to-get-one-back-out
    m = np.zeros(__numel__(ulptr)).flatten().astype(py_real_str)
    p_m = ffi.cast(c_real_str+" *", ffi.from_buffer(m))
    lib.to_numpy(ulptr, p_m)
    m = m.reshape((__rows__(ulptr), __cols__(ulptr)))
    return m
