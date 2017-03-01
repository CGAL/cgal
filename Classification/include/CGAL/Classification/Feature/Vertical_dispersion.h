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
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_FEATURE_VERTICAL_DISPERSION_H
#define CGAL_CLASSIFICATION_FEATURE_VERTICAL_DISPERSION_H

#include <vector>

#include <CGAL/Classification/Image.h>
#include <CGAL/Classification/Planimetric_grid.h>

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

    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
  */
template <typename Geom_traits, typename PointRange, typename PointMap>
class Vertical_dispersion : public Feature_base
{
  typedef Classification::Image<float> Image_float;
  typedef Classification::Planimetric_grid<Geom_traits, PointRange, PointMap> Grid;
  std::vector<double> vertical_dispersion;
  
public:
  /*!
    \brief Constructs the feature.

    \param input input range.
    \param point_map property map to access the input points.
    \param grid precomputed `Planimetric_grid`.
    \param grid_resolution resolution of the planimetric grid.
    \param radius_neighbors radius of local neighborhoods.
  */
  Vertical_dispersion (const PointRange& input,
                       PointMap point_map,
                       const Grid& grid,
                       const double grid_resolution,
                       double radius_neighbors = -1.)
  {
    this->set_weight(1.);
    if (radius_neighbors < 0.)
      radius_neighbors = 5. * grid_resolution;
    
    Image_float Dispersion(grid.width(), grid.height());
    for (std::size_t j = 0; j < grid.height(); j++)	
      for (std::size_t i = 0; i < grid.width(); i++)
        Dispersion(i,j)=0;
    
    std::size_t square = (std::size_t)(0.5 * radius_neighbors / grid_resolution) + 1;
    typename Geom_traits::Vector_3 verti (0., 0., 1.);
    
    for (std::size_t j = 0; j < grid.height(); j++){	
      for (std::size_t i = 0; i < grid.width(); i++){
						
        if(!(grid.mask(i,j)))
          continue;
        std::vector<double> hori;
            
        std::size_t squareXmin = (i < square ? 0 : i - square);
        std::size_t squareXmax = (std::min) (grid.width()-1, i + square);
        std::size_t squareYmin = (j < square ? 0 : j - square);
        std::size_t squareYmax = (std::min) (grid.height()-1, j + square);

        for(std::size_t k = squareXmin; k <= squareXmax; k++)
          for(std::size_t l = squareYmin; l <= squareYmax; l++)
            if(CGAL::sqrt(pow((double)k-i,2)+pow((double)l-j,2))
               <=(double)0.5*radius_neighbors/grid_resolution
               && (grid.indices(k,l).size()>0))
              for(std::size_t t=0; t<grid.indices(k,l).size();t++)
                {
                  std::size_t ip = grid.indices(k,l)[t];
                  hori.push_back (get(point_map, *(input.begin()+ip)).z());
                }
        if (hori.empty())
          continue;
              
        std::sort (hori.begin(), hori.end());

        std::size_t nb_layers = 1;

        std::vector<bool> occupy (1 + (std::size_t)((hori.back() - hori.front()) / grid_resolution), false);
              
        std::size_t last_index = 0;
        for (std::size_t k = 0; k < hori.size(); ++ k)
          {
            std::size_t index = (std::size_t)((hori[k] - hori.front()) / grid_resolution);
            occupy[index] = true;
            if (index > last_index + 1)
              ++ nb_layers;
            last_index = index;
          }

        std::size_t nb_occ = 0;
        for (std::size_t k = 0; k < occupy.size(); ++ k)
          if (occupy[k])
            ++ nb_occ;
					
        Dispersion(i,j)= 1.f - (nb_occ / (float)(occupy.size()));
			
      }
		
    }
    for (std::size_t i = 0; i < input.size(); i++)
      {
        std::size_t I= grid.x(i);
        std::size_t J= grid.y(i);
        vertical_dispersion.push_back((double)Dispersion(I,J));
      }

    this->compute_mean_max (vertical_dispersion, this->mean, this->max);
  }
  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return vertical_dispersion[pt_index];
  }

  virtual std::string name() { return "vertical_dispersion"; }
  /// \endcond
};

}
  
}

}

#endif // CGAL_CLASSIFICATION_FEATURE_VERTICAL_DISPERSION_H
