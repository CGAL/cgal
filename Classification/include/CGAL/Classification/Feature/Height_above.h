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

#ifndef CGAL_CLASSIFICATION_FEATURE_HEIGHT_ABOVE_H
#define CGAL_CLASSIFICATION_FEATURE_HEIGHT_ABOVE_H

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

    %Feature based on local height distribution This feature computes
    the distance between the maximum height on the local cell of the
    planimetric grid and a point's height.

    Its default name is "height_above".

    \tparam GeomTraits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator` and its value type is the key type of
    `PointMap`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `GeomTraits::Point_3`.

  */
template <typename GeomTraits, typename PointRange, typename PointMap>
class Height_above : public Feature_base
{
  using Image_float = Image<float>;
  using Grid = Planimetric_grid<GeomTraits, PointRange, PointMap>;

  const PointRange& input;
  PointMap point_map;
  const Grid& grid;
  Image_float dtm;
  std::vector<float> values;

public:
  /*!
    \brief constructs the feature.

    \param input point range.
    \param point_map property map to access the input points.
    \param grid precomputed `Planimetric_grid`.
  */
  Height_above (const PointRange& input,
                PointMap point_map,
                const Grid& grid)
    : input(input), point_map(point_map), grid(grid)
  {
    this->set_name ("height_above");

    dtm = Image_float(grid.width(),grid.height());

    for (std::size_t j = 0; j < grid.height(); ++ j)
      for (std::size_t i = 0; i < grid.width(); ++ i)
        if (grid.has_points(i,j))
        {
          float z_max = -(std::numeric_limits<float>::max)();

          typename Grid::iterator end = grid.indices_end(i,j);
          for (typename Grid::iterator it = grid.indices_begin(i,j); it != end; ++ it)
          {
            float z = float(get(point_map, *(input.begin()+(*it))).z());
            z_max = ((std::max)(z_max, z));
          }

          dtm(i,j) = z_max;
        }

    if (grid.width() * grid.height() > input.size())
    {
      values.resize (input.size(), 0.f);
      for (std::size_t i = 0; i < input.size(); ++ i)
      {
        std::size_t I = grid.x(i);
        std::size_t J = grid.y(i);
        values[i] = float(dtm(I,J) - get (point_map, *(input.begin() + i)).z());
      }
      dtm.free();
    }

  }

  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t pt_index)
  {
    if (values.empty())
    {
      std::size_t I = grid.x(pt_index);
      std::size_t J = grid.y(pt_index);
      return dtm(I,J) - float(get (point_map, *(input.begin() + pt_index)).z());
    }

    return values[pt_index];
  }

  /// \endcond
};

} // namespace Feature

} // namespace Classification


} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_HEIGHT_ABOVE_H
