// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Florent Lafarge, Simon Giraudot

#ifndef CGAL_CLASSIFICATION_FEATURE_ELEVATION_H
#define CGAL_CLASSIFICATION_FEATURE_ELEVATION_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <CGAL/Classification/Feature_base.h>
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
  typedef typename GeomTraits::Iso_cuboid_3 Iso_cuboid_3;

  typedef Image<float> Image_float;
  typedef Planimetric_grid<GeomTraits, PointRange, PointMap> Grid;

#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
  std::vector<float> elevation_feature;
#else
  const PointRange& input;
  PointMap point_map;
  const Grid& grid;
  Image_float dtm;
#endif
  
public:
  /*!
    \brief Constructs the feature.

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
#ifndef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    : input(input), point_map(point_map), grid(grid)
#endif
  {
    this->set_name ("elevation");
    if (radius_dtm < 0.)
      radius_dtm = 100.f * grid.resolution();

    //DEM
    Image_float dem(grid.width(),grid.height());

    for (std::size_t j = 0; j < grid.height(); ++ j)
      for (std::size_t i = 0; i < grid.width(); ++ i)
        if (grid.has_points(i,j))
        {
          float mean = 0.;
          std::size_t nb = 0;
          typename Grid::iterator end = grid.indices_end(i,j);
          for (typename Grid::iterator it = grid.indices_begin(i,j); it != end; ++ it)
          {
            mean += float(get(point_map, *(input.begin()+(*it))).z());
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

#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    Image_float dtm(grid.width(),grid.height());
#else
    dtm = Image_float(grid.width(),grid.height());
#endif
    
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
          dtm(i,j) = z[z.size() / 10];
        }
    dtm_x.free();

#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    elevation_feature.reserve(input.size());
    for (std::size_t i = 0; i < input.size(); i++){
      std::size_t I = grid.x(i);
      std::size_t J = grid.y(i);
      elevation_feature.push_back ((float)(get(point_map, *(input.begin()+i)).z()-dtm(I,J)));
    }
#endif
    
  }

  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t pt_index)
  {
#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    return elevation_feature[pt_index];
#else
    std::size_t I = grid.x(pt_index);
    std::size_t J = grid.y(pt_index);
    return ((float)(get(point_map, *(input.begin()+pt_index)).z()-dtm(I,J)));
#endif
  }

  /// \endcond
};

} // namespace Feature

} // namespace Classification


} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_ELEVATION_H
