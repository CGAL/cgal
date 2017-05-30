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
//
// Author(s)     : Simon Giraudot, Florent Lafarge

#ifndef CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H
#define CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H

#include <vector>

#include <CGAL/Classification/Image.h>

namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassificationDataStructures

    \brief Class that precomputes a 2D planimetric grid.

    The grid is composed of squared cells with a user-defined size,
    each cell containing the list of indices of the points whose
    projection along the Z-axis lies within this cell. The mapping
    from each point to the cell it lies in is also stored.

    \tparam GeomTraits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator` and its value type is the key type of
    `PointMap`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `GeomTraits::Point_3`.
  */

template <typename GeomTraits, typename PointRange, typename PointMap>
class Planimetric_grid
{
public:
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Iso_cuboid_3 Iso_cuboid_3;

private:
  typedef Image<std::vector<std::size_t> > Image_indices;
  typedef Image<bool> Image_bool;

  Image_indices m_grid;
  float m_resolution;
  std::vector<std::size_t> m_x;
  std::vector<std::size_t> m_y;
  Planimetric_grid* m_lower_scale;
  
public:

  /// \cond SKIP_IN_MANUAL
  Planimetric_grid () { }
  /// \endcond

  /*!
    \brief Constructs a planimetric grid based on the input range.

    \param input point range.
    \param point_map property map to access the input points.
    \param bbox bounding box of the input range.
    \param grid_resolution resolution of the planimetric grid.
  */
  Planimetric_grid (const PointRange& input,
                    PointMap point_map,
                    const Iso_cuboid_3& bbox,
                    float grid_resolution)
    : m_resolution (grid_resolution), m_lower_scale(NULL)
  {
    std::size_t width = (std::size_t)((bbox.xmax() - bbox.xmin()) / grid_resolution) + 1;
    std::size_t height = (std::size_t)((bbox.ymax() - bbox.ymin()) / grid_resolution) + 1;

    m_grid = Image_indices (width, height);

    for (std::size_t i = 0; i < input.size(); ++ i)
    {
      const Point_3& p = get(point_map, *(input.begin()+i));
      m_x.push_back ((std::size_t)((p.x() - bbox.xmin()) / grid_resolution));
      m_y.push_back ((std::size_t)((p.y() - bbox.ymin()) / grid_resolution));
      m_grid(m_x.back(), m_y.back()).push_back (i);
    }
//    std::cerr << "Grid size = " << width << " " << height << std::endl;
  }

  /// \cond SKIP_IN_MANUAL
  Planimetric_grid (Planimetric_grid* lower_scale)
    : m_resolution (lower_scale->resolution() * 2), m_lower_scale (lower_scale)
  {
//    std::cerr << "Grid size = " << width() << " " << height() << std::endl;
  }
  /// \endcond


  /*!
    \brief Returns the resolution of the grid.
  */
  float resolution() const
  {
    return m_resolution;
  }
  
  /*!
    \brief Returns the number of cells along the X-axis.
  */
  std::size_t width() const
  {
    if (m_lower_scale == NULL)
      return m_grid.width();
    else
      return (m_lower_scale->width() + 1) / 2;
  }
  /*!
    \brief Returns the number of cells along the Y-axis.
  */
  std::size_t height() const
  {
    if (m_lower_scale == NULL)
      return m_grid.height();
    else
      return (m_lower_scale->height() + 1) / 2;
  }

  /*!
    \brief Stores the indices of the points lying in the cell at
    position `(x,y)` in `output`.
  */
  template <typename OutputIterator>
  void indices(std::size_t x, std::size_t y, OutputIterator output) const
  {
    if (m_lower_scale == NULL)
    {
      if (x >= m_grid.width() || y >= m_grid.height())
        return;
      std::copy (m_grid(x,y).begin(), m_grid(x,y).end(), output);
    }
    else
    {
      m_lower_scale->indices(x*2, y*2, output);
      m_lower_scale->indices(x*2, y*2 + 1, output);
      m_lower_scale->indices(x*2 + 1, y*2 + 1, output);
      m_lower_scale->indices(x*2 + 1, y*2, output);
    }
  }
  
  /*!
    \brief Returns `false` if the cell at position `(x,y)` is empty, `true` otherwise.
  */
  bool has_points(std::size_t x, std::size_t y) const
  {
    if (m_lower_scale == NULL)
    {
      if (x >= m_grid.width() || y >= m_grid.height())
        return false;
      return (!(m_grid(x,y).empty()));
    }
    else
      return (m_lower_scale->has_points(x*2, y*2)
              || m_lower_scale->has_points(x*2, y*2 + 1)
              || m_lower_scale->has_points(x*2 + 1, y*2 + 1)
              || m_lower_scale->has_points(x*2 + 1, y*2));
  }

  /*!
    \brief Returns the `x` grid coordinate of the point at position `index`.
  */
  std::size_t x(std::size_t index) const
  {
    if (m_lower_scale == NULL)
      return m_x[index];
    else
      return m_lower_scale->x(index) / 2;
  }
  /*!
    \brief Returns the `y` grid coordinate of the point at position `index`.
  */
  std::size_t y(std::size_t index) const
  {
    if (m_lower_scale == NULL)
      return m_y[index];
    else
      return m_lower_scale->y(index) / 2;
  }
};
  

}
  
}


#endif // CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H
