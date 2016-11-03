// Copyright (c) 2016  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_CLASSIFICATION_ATTRIBUTE_ELEVATION_H
#define CGAL_CLASSIFICATION_ATTRIBUTE_ELEVATION_H

#include <vector>

#include <CGAL/Classification/Attribute.h>
#include <CGAL/Classification/Image.h>
#include <CGAL/Classification/Planimetric_grid.h>

namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassification

    \brief Attribute based on local elevation.

    The local position of the ground can be computed for urban
    scenes. This attribute computes the distance to the local
    estimation of the ground.

    It is useful to discriminate the ground from horizontal roofs.

    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.

  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap>
class Attribute_elevation : public Attribute
{
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;

  typedef Image<float> Image_float;
  typedef Planimetric_grid<Kernel, RandomAccessIterator, PointMap> Grid;
   
  std::vector<double> elevation_attribute;
  
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points
    \param grid Precomputed `Planimetric_grid`
    \param grid_resolution Resolution of the planimetric grid
    \param radius_dtm Radius for digital terrain modeling (must be bigger than the size of a building)

  */
  Attribute_elevation (RandomAccessIterator begin,
                       RandomAccessIterator end,
                       PointMap point_map,
                       const Grid& grid,
                       const double grid_resolution,
                       double radius_dtm = -1.)
  {
    this->weight = 1.;
    if (radius_dtm < 0.)
      radius_dtm = 100. * grid_resolution;

    //DEM
    Image_float dem(grid.width(),grid.height());

    for (std::size_t j = 0; j < grid.height(); ++ j)
      for (std::size_t i = 0; i < grid.width(); ++ i)
        {
          if (grid.indices(i,j).empty())
            continue;
          double mean = 0.;
          for (std::size_t k = 0; k < grid.indices(i,j).size(); ++ k)
            mean += get(point_map, begin[grid.indices(i,j)[k]]).z();
          mean /= grid.indices(i,j).size();
          dem(i,j) = mean;
        }

    std::size_t square = (std::size_t)(0.5 * radius_dtm / grid_resolution) + 1;
    
    Image_float dtm_x(grid.width(),grid.height());
    
    for (std::size_t j = 0; j < grid.height(); ++ j)
      {
        for (std::size_t i = 0; i < grid.width(); ++ i)
          {
            std::size_t squareXmin = (i < square ? 0 : i - square);
            std::size_t squareXmax = (std::min)(grid.width() - 1, i + square);

            std::vector<double> z;
            z.reserve(squareXmax - squareXmin +1 );
            for(std::size_t k = squareXmin; k <= squareXmax; k++)
              if (dem(k,j) != 0.)
                z.push_back (dem(k,j));
            if (z.empty())
              continue;
            std::nth_element (z.begin(), z.begin() + (z.size() / 10), z.end());
            dtm_x(i,j) = z[z.size() / 10];
          }
      }
    dem.free();

    Image_float dtm(grid.width(),grid.height());
    
    for (std::size_t i = 0; i < grid.width(); ++ i)
      {
        for (std::size_t j = 0; j < grid.height(); ++ j)
          {
            std::size_t squareYmin = (j < square ? 0 : j - square);
            std::size_t squareYmax = (std::min)(grid.height() - 1, j + square);
            std::vector<double> z;
            z.reserve(squareYmax - squareYmin +1 );
            for(std::size_t l = squareYmin; l <= squareYmax; l++)
              if (dtm_x(i,l) != 0.)
                z.push_back (dtm_x(i,l));
            if (z.empty())
              continue;
            std::nth_element (z.begin(), z.begin() + (z.size() / 10), z.end());
            dtm(i,j) = z[z.size() / 10];
          }
      }
    dtm_x.free();

    elevation_attribute.reserve(end - begin);
    for (std::size_t i = 0; i < (std::size_t)(end - begin); i++){
      std::size_t I = grid.x(i);
      std::size_t J = grid.y(i);
      elevation_attribute.push_back ((double)(get(point_map, begin[i]).z()-dtm(I,J)));
    }

    this->compute_mean_max (elevation_attribute, this->mean, this->max);
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return elevation_attribute[pt_index];
  }
  
  virtual std::string id() { return "elevation"; }
  /// \endcond
};

} // namespace Classification


} // namespace CGAL

#endif // CGAL_CLASSIFICATION_ATTRIBUTE_ELEVATION_H
