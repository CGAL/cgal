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
// Author(s)     : Simon Giraudot, Florent Lafarge

#ifndef CGAL_CLASSIFICATION_FEATURE_ECHO_SCATTER_H
#define CGAL_CLASSIFICATION_FEATURE_ECHO_SCATTER_H

#include <CGAL/license/Classification.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/number_utils.h>
#include <vector>
#include <cmath>


namespace CGAL {

namespace Classification {

namespace Feature {
  
  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on echo scatter. The number of returns (echo
    number) is a useful information provided by most LIDAR sensors. It
    can help to identify trees.

    Its default name is "echo_scatter".

    \tparam GeomTraits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator` and its value type is the key type of
    `PointMap`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `GeomTraits::Point_3`.
    \tparam EchoMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `std::size_t`.
  */
template <typename GeomTraits, typename PointRange, typename PointMap, typename EchoMap>
class Echo_scatter : public Feature_base
{
public:
  typedef Classification::Planimetric_grid<GeomTraits, PointRange, PointMap> Grid;
private:  
  typedef Classification::Image<float> Image_float;

  std::vector<float> echo_scatter;
  
public:
  /*!
    \brief Constructs the feature.

    \param input point range.
    \param echo_map property map to access the echo values of the input points.
    \param grid precomputed `Planimetric_grid`.
    \param radius_neighbors radius of local neighborhoods.
  */
  Echo_scatter (const PointRange& input,
                EchoMap echo_map,
                const Grid& grid,
                float radius_neighbors = 1.)
  {
    this->set_name ("echo_scatter");
    Image_float Scatter(grid.width(), grid.height());
    for (std::size_t j = 0; j < grid.height(); j++)
      for (std::size_t i = 0; i < grid.width(); i++)
        if (grid.has_points(i,j))
          Scatter(i,j)=0;

    std::size_t square = (std::size_t)(0.5 * radius_neighbors / grid.resolution()) + 1;

    for (std::size_t j = 0; j < grid.height(); j++)
      for (std::size_t i = 0; i < grid.width(); i++)
        if(grid.has_points(i,j))
        {

          std::size_t squareXmin = (i < square ? 0 : i - square);
          std::size_t squareXmax = (std::min) (grid.width()-1, i + square);
          std::size_t squareYmin = (j < square ? 0 : j - square);
          std::size_t squareYmax = (std::min) (grid.height()-1, j + square);
			
          std::size_t NB_echo_sup=0;
          std::size_t NB_echo_total=0;

          for(std::size_t k = squareXmin; k <= squareXmax; k++){
            for(std::size_t l = squareYmin; l <= squareYmax; l++){
									
              if(CGAL::sqrt(pow((float)k-i,2)+pow((float)l-j,2))<=(float)0.5*radius_neighbors/grid.resolution())
              {
                typename Grid::iterator end = grid.indices_end(k,l);
                std::size_t nb = 0;
                for (typename Grid::iterator it = grid.indices_begin(k,l); it != end; ++ it)
                {
                  ++ nb;
                  if(get(echo_map, *(input.begin()+(*it))) > 1)
                    NB_echo_sup++;
                }

                NB_echo_total=NB_echo_total+nb;

              }
						
            }
					
          }
					
          Scatter(i,j)=(float)NB_echo_sup/NB_echo_total;
        }
  
    for(std::size_t i = 0; i < input.size(); i++){
      std::size_t I= grid.x(i);
      std::size_t J= grid.y(i);
      echo_scatter.push_back((float)Scatter(I,J));
    }
  }

  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t pt_index)
  {
    return echo_scatter[pt_index];
  }
  /// \endcond
};

} // namespace Feature

} // namespace Classification
  
} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_ECHO_SCATTER_H
