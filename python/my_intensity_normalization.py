import numpy as np
from scipy.stats import norm
from scipy.ndimage import gaussian_filter


def my_intensity_normalization(struct_img, scaling_param, verbose=False):

    '''
    Mode 1:  scaling_param = [0]
    Mode 2:  scaling_param = [lower std range, upper std range]
    Mode 3:  scaling_param = [lower std range, upper std range, lower abs intensity, higher abs intensity]
    '''
    assert len(scaling_param) > 0

    if len(scaling_param) == 1:
        if scaling_param[0] < 1:
            if verbose: print('intensity normalization: using min-max normalization with NO absolute intensity upper bound')
        else:
            if verbose: print(f'intensity normalization: using min-max normalization with absolute intensity upper bound {scaling_param[0]}')
            struct_img[struct_img > scaling_param[0]] = struct_img.min()
        strech_min = struct_img.min()
        strech_max = struct_img.max()
        struct_img = (struct_img - strech_min + 1e-8)/(strech_max - strech_min + 1e-8)
    elif len(scaling_param) == 2:
        if verbose:  print(f'intensity normalization: normalize into [mean - {scaling_param[0]} x std, mean + {scaling_param[1]} x std] ')
        m, s = norm.fit(struct_img.flat)
        strech_min = max(m - scaling_param[0] * s, struct_img.min())
        strech_max = min(m + scaling_param[1] * s, struct_img.max())
        struct_img[struct_img > strech_max] = strech_max
        struct_img[struct_img < strech_min] = strech_min
        struct_img = (struct_img - strech_min + 1e-8)/(strech_max - strech_min + 1e-8)
    elif len(scaling_param) == 4:
        img_valid = struct_img[np.logical_and(struct_img > scaling_param[2], struct_img < scaling_param[3])]
        m, s = norm.fit(img_valid.flat)
        strech_min = max(scaling_param[2] - scaling_param[0] * s, struct_img.min())
        strech_max = min(scaling_param[3] + scaling_param[1] * s, struct_img.max())
        struct_img[struct_img > strech_max] = strech_max
        struct_img[struct_img < strech_min] = strech_min
        struct_img = (struct_img - strech_min + 1e-8)/(strech_max - strech_min + 1e-8)

    if verbose:  print('intensity normalization completes')
    return struct_img
