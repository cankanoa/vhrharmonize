"""Explicit lossless TIFF compression selection for native raster plugins."""
def validate_compression(compression):
    if not isinstance(compression, bool):
        raise ValueError("compression must be true or false")
    return compression


def raster_profile(profile, compression=True):
    """Return a copy; source TIFF compression must not leak into uncompressed output."""
    validate_compression(compression)
    result = dict(profile)
    if compression:
        result['compress'] = result.get('compress') or 'deflate'
    else:
        for key in ('compress', 'predictor', 'zlevel', 'zstd_level', 'jpeg_quality', 'level'):
            result.pop(key, None)
        result['compress'] = 'none'
    return result
