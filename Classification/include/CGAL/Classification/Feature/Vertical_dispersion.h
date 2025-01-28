// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_FEATURE_VERTICAL_DISPERSION_H
#define CGAL_CLASSIFICATION_FEATURE_VERTICAL_DISPERSION_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <CGAL/number_utils.h>
#include <CGAL/Classification/Image.h>
#include <CGAL/Classification/compressed_float.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <boost/algorithm/minmax_element.hpp>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/int.h>
#include <CGAL/float.h>
#include <boost/tuple/tuple.hpp>

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on local vertical dispersion of points. Urban
    scenes can often be decomposed as a set of 2D regions with
    different heights. While these heights are usually piecewise
    constant or piecewise linear, on some specific parts of the scene
    such as vegetation, they can become extremely unstable. This
    feature quantifies the vertical dispersion of the points on a
    local Z-cylinder around the points.

    Its default name is "vertical_dispersion".

    \tparam GeomTraits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator` and its value type is the key type of
    `PointMap`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `GeomTraits::Point_3`.
  */
template <typename GeomTraits, typename PointRange, typename PointMap>
class Vertical_dispersion : public Feature_base
{
  using Image_cfloat = Classification::Image<compressed_float>;
  using Grid = Classification::Planimetric_grid<GeomTraits, PointRange, PointMap>;

  const Grid& grid;
  Image_cfloat Dispersion;
  std::vector<compressed_float> values;

public:
  /*!
    \brief constructs the feature.

    \param input point range.
    \param point_map property map to access the input points.
    \param grid precomputed `Planimetric_grid`.
    \param radius_neighbors radius of local neighborhoods.
  */
  Vertical_dispersion (const PointRange& input,
                       PointMap point_map,
                       const Grid& grid,
                       float radius_neighbors = -1.)
    : grid (grid)
  {
    this->set_name ("vertical_dispersion");
    if (radius_neighbors < 0.)
      radius_neighbors = 5.f * grid.resolution();

    if (grid.width() * grid.height() > input.size())
      values.resize(input.size(), compressed_float(0));
    else
    {
      Dispersion = Image_cfloat(grid.width(), grid.height());
      for (std::size_t j = 0; j < grid.height(); j++)
        for (std::size_t i = 0; i < grid.width(); i++)
          if (grid.has_points(i,j))
            Dispersion(i,j) = compressed_float(0);
    }

    std::size_t square = (std::size_t)(0.5 * radius_neighbors / grid.resolution()) + 1;
    typename GeomTraits::Vector_3 verti (0., 0., 1.);

    std::vector<float> hori;

    for (std::size_t j = 0; j < grid.height(); j++){
      for (std::size_t i = 0; i < grid.width(); i++){

        if(!(grid.has_points(i,j)))
          continue;
        hori.clear();

        std::size_t squareXmin = (i < square ? 0 : i - square);
        std::size_t squareXmax = (std::min) (grid.width()-1, i + square);
        std::size_t squareYmin = (j < square ? 0 : j - square);
        std::size_t squareYmax = (std::min) (grid.height()-1, j + square);

        float bound = (float)0.5*radius_neighbors/grid.resolution();
        bound = CGAL::square(bound);
        for(std::size_t k = squareXmin; k <= squareXmax; k++)
          for(std::size_t l = squareYmin; l <= squareYmax; l++)
          {
            if(CGAL::square((float)(k-i))+ CGAL::square((float)(l-j))
               <= bound)
            {
              for (typename Grid::iterator it = grid.indices_begin(k,l); it != grid.indices_end(k,l); ++ it)
                hori.push_back (float(get(point_map, *(input.begin()+(*it))).z()));
            }
          }

        if (hori.empty())
          continue;

        std::vector<float>::iterator min_it, max_it;
        boost::tie(min_it, max_it)
          = boost::minmax_element (hori.begin(), hori.end());

        std::vector<bool> occupy (1 + (std::size_t)((*max_it - *min_it) / grid.resolution()), false);

        for (std::size_t k = 0; k < hori.size(); ++ k)
        {
          std::size_t index = (std::size_t)((hori[k] - *min_it) / grid.resolution());
          occupy[index] = true;
        }

        std::size_t nb_occ = 0;
        for (std::size_t k = 0; k < occupy.size(); ++ k)
          if (occupy[k])
            ++ nb_occ;

        compressed_float v = compress_float (1.f - (nb_occ / (float)(occupy.size())));
        if (values.empty())
          Dispersion(i,j) = v;
        else
        {
          typename Grid::iterator end = grid.indices_end(i,j);
          for (typename Grid::iterator it = grid.indices_begin(i,j); it != end; ++ it)
            values[*it] = v;
        }
      }

    }
  }
  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t pt_index)
  {
    if (values.empty())
    {
      std::size_t I = grid.x(pt_index);
      std::size_t J = grid.y(pt_index);
      return decompress_float (Dispersion(I,J));
    }
    return decompress_float (values[pt_index]);
  }
  /// \endcond
};

}

}

}

#endif // CGAL_CLASSIFICATION_FEATURE_VERTICAL_DISPERSION_H
