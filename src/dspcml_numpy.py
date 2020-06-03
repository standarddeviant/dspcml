import numpy as np
from _dspcml_numpy import ffi, lib

py_real_str = "float32"
c_real_str = "float"


def print_matrix(ulptr):
    strarg = ffi.new("char[4096]")
    lib.print_matrix(ulptr, strarg)
    print(ffi.string(strarg).decode())


def maybe_print_errno():
    strarg = ffi.new("char[4096]")
    lib.maybe_print_errno(strarg)
    print(ffi.string(strarg).decode())


def zeros(rows, cols):
    # return DSPCML MATRIX pointer as an unsigned long for simplicity
    out = lib.zeros(rows, cols)
    maybe_print_errno()
    return out


def arange(start, stop, step):
    return lib.arange(start, stop, step)


def tile(ulptr, row_copies, col_copies):
    return lib.tile(ulptr, row_copies, col_copies)


def new_sos(nb_bq, nb_ch, coeffs):
    arr = coeffs.flatten(order='C').astype(py_real_str)
    p_arr = ffi.cast(c_real_str+" *", ffi.from_buffer(arr))
    return lib.new_sos(nb_bq, nb_ch, p_arr)


def cpy_submat(ulptr_dst, ulptr_src,
               dst_row, dst_col,
               src_row, src_col,
               nb_rows, nb_cols):
    lib.cpy_submat(ulptr_dst, ulptr_src,
                   dst_row, dst_col,
                   src_row, src_col,
                   nb_rows, nb_cols)


def __numel__(ulptr):
    # return number of elements in DSPCML MATRIX
    return lib.numel(ulptr)


def __rows__(ulptr):
    return lib.rows(ulptr)


def __cols__(ulptr):
    return lib.cols(ulptr)


def __sin__(ulptr, ulptr_opt):
    return lib.__sin__(ulptr, ulptr_opt)


def __cos__(ulptr, ulptr_opt):
    return lib.__cos__(ulptr, ulptr_opt)


def reshape(ulptr, new_rows, new_cols):
    return lib.reshape(ulptr, new_rows, new_cols)


def to_numpy(ulptr):
    # Use numpy/ctypes trick to make valid C pointer from numpy array
    # https://stackoverflow.com/questions/16276268/how-to-pass-a-numpy-array-into-a-cffi-function-and-how-to-get-one-back-out
    m = np.zeros(__numel__(ulptr)).flatten().astype(py_real_str)
    p_m = ffi.cast(c_real_str+" *", ffi.from_buffer(m))
    lib.to_numpy(ulptr, p_m)
    m = m.reshape((__rows__(ulptr), __cols__(ulptr)))
    return m


def to_dspcml(arr):
    if len(arr.shape) > 2:
        pass # TODO throw error

    arr = arr.astype(py_real_str)
    p_arr = ffi.cast(c_real_str+" *", ffi.from_buffer(arr))

    rows = arr.shape[0]
    cols = 1
    if 2 == len(arr.shape):
        cols = arr.shape[1]

    ulptr = lib.to_dspcml(p_arr, rows, cols)

    return ulptr
