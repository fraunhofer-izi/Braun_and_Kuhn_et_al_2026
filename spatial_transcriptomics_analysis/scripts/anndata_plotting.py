# Copyright 2026 Fraunhofer-Gesellschaft zur Förderung der angewandten
# Forschung e.V.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.


import scanpy as sc
import logging

logger = logging.getLogger(__name__)


def crop_region_by_fraction(adata, xfrac=(0, 1), yfrac=(0, 1)):
    """
    Crop bounding box based on fractions of the x and y ranges.

    Parameters
    ----------
    adata : AnnData
        Object with .obsm["spatial"] coordinates.
    xfrac : tuple (min_frac, max_frac)
        Fraction of x-range to keep. (0,1) = full, (0,0.5) = left half, (0.5,1) = right half.
    yfrac : tuple (min_frac, max_frac)
        Fraction of y-range to keep. (0,1) = full, (0,0.5) = bottom half, (0.5,1) = top half.

    Returns
    -------
    crop_box : tuple (xmin, xmax, ymin, ymax)
    """
    coords = adata.obsm["spatial"]
    xmin, ymin = coords.min(axis=0)
    xmax, ymax = coords.max(axis=0)

    xrange = xmax - xmin
    yrange = ymax - ymin

    xmin_new = xmin + xrange * xfrac[0]
    xmax_new = xmin + xrange * xfrac[1]
    ymin_new = ymin + yrange * yfrac[0]
    ymax_new = ymin + yrange * yfrac[1]

    return (xmin_new, xmax_new, ymin_new, ymax_new)
