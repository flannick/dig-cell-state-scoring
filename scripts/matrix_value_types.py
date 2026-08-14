import numpy as np

VALUE_TYPES = {'raw_counts', 'linear_cp10k', 'log1p_cp10k', 'linear_normalized', 'log1p_normalized'}


def require_nonnegative(has_negative, value_type):
    if has_negative:
        raise SystemExit(f'Matrix has negative values but --value-type={value_type}')


def linearize_expression_matrix(matrix, value_type):
    mat = matrix.astype(float).tocsr(copy=True)
    if value_type == 'raw_counts':
        totals = np.asarray(mat.sum(axis=1)).ravel()
        return mat.multiply(10000.0 / totals[:, None]).tocsr()
    if value_type in ('linear_cp10k', 'linear_normalized'):
        return mat
    mat.data = np.expm1(mat.data)
    mat.data[mat.data < 0] = 0.0
    return mat.tocsr()
