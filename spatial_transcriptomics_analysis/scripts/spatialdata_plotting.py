# Copyright 2026 Fraunhofer-Gesellschaft zur Förderung der angewandten
# Forschung e.V.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.

from spatialdata import bounding_box_query
import logging

logger = logging.getLogger(__name__)


def crop_sdata_bounding_box(
    sdata,
    min_coordinate=None,
    max_coordinate=None,
    target_coordinate_system=None,
):
    """
    Crop a SpatialData object using an axis-aligned bounding box.

    This is a thin wrapper around ``spatialdata.bounding_box_query`` that
    standardizes axis usage for 2D spatial data.

    Parameters
    ----------
    sdata : SpatialData
        SpatialData object to crop.
    min_coordinate : tuple[float, float] or None, optional
        Lower (x, y) corner of the bounding box.
        If None, no lower bound is applied.
    max_coordinate : tuple[float, float] or None, optional
        Upper (x, y) corner of the bounding box.
        If None, no upper bound is applied.
    target_coordinate_system : str or None, optional
        Coordinate system in which the crop is performed.

    Returns
    -------
    SpatialData
        Cropped SpatialData object.
    """
    return bounding_box_query(
        sdata,
        min_coordinate=min_coordinate,
        max_coordinate=max_coordinate,
        axes=("x", "y"),
        target_coordinate_system=target_coordinate_system,
    )


def crop_region_by_fraction_from_shape_layer(
    sdata, cs_name, shape_layer, xfrac=(0, 1), yfrac=(0, 1)
):
    """
    Crop a SpatialData object using fractional bounds of a shape layer.

    The bounding box of the specified shape layer is computed, and a
    sub-region is retained based on fractional ranges along x and y.

    Parameters
    ----------
    sdata : SpatialData
        SpatialData object to crop.
    cs_name : str
        Coordinate system in which the crop is performed.
    shape_layer : str
        Name of the shape layer whose bounding box is used
        (e.g. ``"cell_segmentation"``).
    xfrac : tuple[float, float], optional
        Fractional range (min, max) of the x-axis to keep.
        Values must lie in [0, 1].
    yfrac : tuple[float, float], optional
        Fractional range (min, max) of the y-axis to keep.
        Values must lie in [0, 1].

    Returns
    -------
    SpatialData
        Cropped SpatialData object.

    Notes
    -----
    This function is useful for selecting central regions or sub-tiles
    of large Visium HD samples in a resolution-independent way.
    """
    # Get bounding box from the shape layer
    min_x, min_y, max_x, max_y = sdata[shape_layer].geometry.total_bounds

    # Compute fractional ranges
    x0 = min_x + (max_x - min_x) * xfrac[0]
    x1 = min_x + (max_x - min_x) * xfrac[1]
    y0 = min_y + (max_y - min_y) * yfrac[0]
    y1 = min_y + (max_y - min_y) * yfrac[1]

    # Crop the spatialdata object
    return crop0(
        sdata,
        min_coordinate=[x0, y0],
        max_coordinate=[x1, y1],
        target_coordinate_system=cs_name,
    )


def crop_sdata_to_shape_bounds(sdata, shape_layer, coordinate_system):
    """
    Crop a SpatialData object to the bounding box of a shape layer.

    The total bounds of the specified shape layer are computed and used
    to crop all spatial elements in the given coordinate system.

    Parameters
    ----------
    sdata : SpatialData
        SpatialData object to crop.
    shape_layer : str
        Name of the shape layer whose bounding box defines the crop
        (e.g. ``"cell_segmentation"``).
    coordinate_system : str
        Coordinate system in which the crop is performed.

    Returns
    -------
    SpatialData
        Cropped SpatialData object.
    """
    # Get bounding box from the specified shape layer
    min_x, min_y, max_x, max_y = sdata[shape_layer].geometry.total_bounds

    # Crop sdata using bounding box
    sdata_cropped = sdata.query.bounding_box(
        min_coordinate=[min_x, min_y],
        max_coordinate=[max_x, max_y],
        axes=("x", "y"),
        target_coordinate_system=coordinate_system,
    )

    return sdata_cropped


def crop0(x, min_coordinate=None, max_coordinate=None, target_coordinate_system=None):
    return bounding_box_query(
        x,
        min_coordinate=min_coordinate,
        max_coordinate=max_coordinate,
        axes=("x", "y"),
        target_coordinate_system=target_coordinate_system,
    )
