// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Florent Lafarge, Simon Giraudot

#ifndef CGAL_CLASSIFICATION_FEATURE_ELEVATION_H
#define CGAL_CLASSIFICATION_FEATURE_ELEVATION_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/compressed_float.h>
#include <CGAL/Classification/Image.h>
#include <CGAL/Classification/Planimetric_grid.h>

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on local elevation. The local position of the
    ground can be computed for urban scenes. This feature computes
    the distance to the local estimation of the ground. It is useful
    to discriminate the ground from horizontal roofs.

    Its default name is "elevation".

    \tparam GeomTraits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator` and its value type is the key type of
    `PointMap`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `GeomTraits::Point_3`.

  */
template <typename GeomTraits, typename PointRange, typename PointMap>
class Elevation : public Feature_base
{
  using Image_float = Image<float>;
  using Image_cfloat = Image<compressed_float>;
  using Grid = Planimetric_grid<GeomTraits, PointRange, PointMap>;

  const PointRange& input;
  PointMap point_map;
  const Grid& grid;
  Image_cfloat dtm;
  std::vector<compressed_float> values;
  float z_max;
  float z_min;

public:
  /*!
    \brief constructs the feature.

    \param input point range.
    \param point_map property map to access the input points.
    \param grid precomputed `Planimetric_grid`.
    \param radius_dtm radius for digital terrain modeling (should be
    larger than the width and length of the largest building).
  */
  Elevation (const PointRange& input,
             PointMap point_map,
             const Grid& grid,
             float radius_dtm = -1.)
    : input(input), point_map(point_map), grid(grid)
  {
    this->set_name ("elevation");
    if (radius_dtm < 0.)
      radius_dtm = 10.f * grid.resolution();

    //DEM
    Image_float dem(grid.width(),grid.height());

    z_max = 0.f;
    z_min = (std::numeric_limits<float>::max)();

    for (std::size_t j = 0; j < grid.height(); ++ j)
      for (std::size_t i = 0; i < grid.width(); ++ i)
        if (grid.has_points(i,j))
        {
          float mean = 0.;
          std::size_t nb = 0;
          typename Grid::iterator end = grid.indices_end(i,j);
          for (typename Grid::iterator it = grid.indices_begin(i,j); it != end; ++ it)
          {
            float z = float(get(point_map, *(input.begin()+(*it))).z());
            z_min = ((std::min)(z_min, z));
            z_max = ((std::max)(z_max, z));
            mean += z;
            ++ nb;
          }
          if (nb == 0)
            continue;
          mean /= nb;
          dem(i,j) = mean;
        }

    std::size_t square = (std::size_t)(0.5 * radius_dtm / grid.resolution()) + 1;

    Image_float dtm_x(grid.width(),grid.height());

    for (std::size_t j = 0; j < grid.height(); ++ j)
      for (std::size_t i = 0; i < grid.width(); ++ i)
        if (grid.has_points(i,j))
        {
          std::size_t squareXmin = (i < square ? 0 : i - square);
          std::size_t squareXmax = (std::min)(grid.width() - 1, i + square);

          std::vector<float> z;
          z.reserve(squareXmax - squareXmin +1 );
          for(std::size_t k = squareXmin; k <= squareXmax; k++)
            if (dem(k,j) != 0.)
              z.push_back (dem(k,j));
          if (z.empty())
            continue;
          std::nth_element (z.begin(), z.begin() + (z.size() / 10), z.end());
          dtm_x(i,j) = z[z.size() / 10];
        }
    dem.free();

    if (grid.width() * grid.height() > input.size())
      values.resize (input.size(), compressed_float(0));
    else
      dtm = Image_cfloat(grid.width(),grid.height());

    for (std::size_t i = 0; i < grid.width(); ++ i)
      for (std::size_t j = 0; j < grid.height(); ++ j)
        if (grid.has_points(i,j))
        {
          std::size_t squareYmin = (j < square ? 0 : j - square);
          std::size_t squareYmax = (std::min)(grid.height() - 1, j + square);
          std::vector<float> z;
          z.reserve(squareYmax - squareYmin +1 );
          for(std::size_t l = squareYmin; l <= squareYmax; l++)
            if (dtm_x(i,l) != 0.)
              z.push_back (dtm_x(i,l));
          if (z.empty())
            continue;
          std::nth_element (z.begin(), z.begin() + (z.size() / 10), z.end());

          compressed_float v = compress_float (z[z.size() / 10], z_min, z_max);
          if (values.empty())
            dtm(i,j) = v;
          else
          {
            typename Grid::iterator end = grid.indices_end(i,j);
            for (typename Grid::iterator it = grid.indices_begin(i,j); it != end; ++ it)
              values[*it] = v;
          }
        }
    dtm_x.free();

  }

  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t pt_index)
  {
    float d = 0.f;
    if (values.empty())
    {
      std::size_t I = grid.x(pt_index);
      std::size_t J = grid.y(pt_index);
      d = decompress_float (dtm(I,J), z_min, z_max);
    }
    else
      d = decompress_float (values[pt_index], z_min, z_max);

    return ((float)(get(point_map, *(input.begin()+pt_index)).z()-d));
  }

  /// \endcond
};

} // namespace Feature

} // namespace Classification


} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_ELEVATION_H
