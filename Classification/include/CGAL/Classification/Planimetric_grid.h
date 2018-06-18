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

#ifndef CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H
#define CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <CGAL/Classification/Image.h>

#include <boost/iterator/iterator_facade.hpp>

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
  std::size_t m_current_scale;

  std::size_t m_width;
  std::size_t m_height;
  std::vector<bool> m_has_points;
  
public:

#ifdef DOXYGEN_RUNNING
  typedef unspecified_type iterator; ///< A forward iterator with value type `std::size_t`.
#else
  class iterator
    : public boost::iterator_facade<iterator, std::size_t, std::forward_iterator_tag>
  {
  public:
    friend class boost::iterator_core_access;
    
    iterator(const Planimetric_grid* lowest_scale,
             std::size_t scale,
             std::size_t large_x,
             std::size_t large_y,
             bool end = false)
      : m_lowest_scale (lowest_scale)
      , m_idx(0)
    {
      std::size_t size = 1;

      while (scale != 0)
      {
        size *= 2;
        -- scale;
      }

      std::size_t xmin = large_x * size;
      m_xmax = (large_x + 1) * size;
      m_ymin = large_y * size;
      m_ymax = (large_y + 1) * size;
      
      m_pos_x = xmin;
      m_pos_y = m_ymin;

      bool found_one = false;
      for (std::size_t x = xmin; x < m_xmax; ++ x)
      {
        for (std::size_t y = m_ymin; y < m_ymax; ++ y)
          if (lowest_scale_has_points(x,y))
          {
            m_pos_x = x;
            m_pos_y = y;
            found_one = true;
            break;
          }
        if(found_one)
          break;
      }

      if (end && found_one)
      {
        m_pos_x = m_xmax;
        m_pos_y = m_ymax;
      }
    }

    void increment()
    {
      ++ m_idx;
      if (m_idx == m_lowest_scale->m_grid(m_pos_x, m_pos_y).size())
      {
        m_idx = 0;
        do
        {
          ++ m_pos_y;
          if (m_pos_y == m_ymax)
          {
            m_pos_y = m_ymin;
            
            ++ m_pos_x;
            if (m_pos_x == m_xmax) // end() reached
            {
              m_pos_y = m_ymax; // put y to max so that this == end()
              break;
            }
              
          }
        }
        while (!(lowest_scale_has_points(m_pos_x, m_pos_y)));
      }
    }

    bool equal (const iterator& other) const
    {
      return (m_pos_x == other.m_pos_x &&
              m_pos_y == other.m_pos_y &&
              m_idx == other.m_idx);
    }

    std::size_t& dereference() const
    {
      return const_cast<std::size_t&>(m_lowest_scale->m_grid(m_pos_x, m_pos_y)[m_idx]);
    }

  private:

    const Planimetric_grid* m_lowest_scale;
    std::size_t m_xmin, m_xmax, m_ymin, m_ymax;

    std::size_t m_size;
    std::size_t m_idx;
    std::size_t m_pos_x;
    std::size_t m_pos_y;

    bool lowest_scale_has_points (std::size_t x, std::size_t y) const
    {
      if (x >= m_lowest_scale->width() || y >= m_lowest_scale->height())
        return false;
      return m_lowest_scale->has_points (x, y);
    }

  };
#endif

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
    : m_resolution (grid_resolution), m_lower_scale(NULL), m_current_scale(0)
  {
    m_width = (std::size_t)((bbox.xmax() - bbox.xmin()) / grid_resolution) + 1;
    m_height = (std::size_t)((bbox.ymax() - bbox.ymin()) / grid_resolution) + 1;

    m_grid = Image_indices (m_width, m_height);

    for (std::size_t i = 0; i < input.size(); ++ i)
    {
      const Point_3& p = get(point_map, *(input.begin()+i));
      m_x.push_back ((std::size_t)((p.x() - bbox.xmin()) / grid_resolution));
      m_y.push_back ((std::size_t)((p.y() - bbox.ymin()) / grid_resolution));
      m_grid(m_x.back(), m_y.back()).push_back (i);
    }
  }

  /// \cond SKIP_IN_MANUAL
  Planimetric_grid (Planimetric_grid* lower_scale)
    : m_resolution (lower_scale->resolution() * 2), m_lower_scale (lower_scale)
  {
    m_current_scale = lower_scale->m_current_scale + 1;

    m_width = (m_lower_scale->width() + 1) / 2;
    m_height = (m_lower_scale->height() + 1) / 2;
    
    m_has_points.reserve(m_width * m_height);
    for (std::size_t x = 0; x < m_width; ++ x)
      for (std::size_t y = 0; y < m_height; ++ y)
      {
        bool has_points = false;

        for (std::size_t i = 0; i <= 1; ++ i)
        {
          std::size_t xi = x*2 + i;
          if (xi >= m_lower_scale->width())
            continue;

          for (std::size_t j = 0; j <= 1; ++ j)
          {
            std::size_t yi = y*2 + j;
            if (yi >= m_lower_scale->height())
              continue;

            if (m_lower_scale->has_points(xi,yi))
            {
              has_points = true;
              break;
            }
          }
          if (has_points)
            break;
        }

        m_has_points.push_back (has_points);
      }
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
    return m_width;
  }
  /*!
    \brief Returns the number of cells along the Y-axis.
  */
  std::size_t height() const
  {
    return m_height;
  }

  /// \cond SKIP_IN_MANUAL
  const Planimetric_grid* lowest_scale() const
  {
    if (m_current_scale == 0)
      return this;

    // else
    return m_lower_scale->lowest_scale();
  }
  /// \endcond

  /*!
    \brief Returns the begin iterator on the indices of the points
    lying in the cell at position `(x,y)`.
  */
  iterator indices_begin(std::size_t x, std::size_t y) const
  {
    CGAL_assertion (x < m_width && y < m_height);

    return iterator (lowest_scale(), m_current_scale, x, y);
  }

  /*!
    \brief Returns the past-the-end iterator on the indices of the points
    lying in the cell at position `(x,y)`.
  */
  iterator indices_end(std::size_t x, std::size_t y) const
  {
    CGAL_assertion (x < m_width && y < m_height);

    return iterator (lowest_scale(), m_current_scale, x, y, true);
  }

  /*!
    \brief Returns `false` if the cell at position `(x,y)` is empty, `true` otherwise.
  */
  bool has_points(std::size_t x, std::size_t y) const
  {
    CGAL_assertion (x < m_width && y < m_height);
    
    if (m_current_scale == 0)
      return (!(m_grid(x,y).empty()));

    // else
    return m_has_points[x * m_height + y];
  }

  /*!
    \brief Returns the `x` grid coordinate of the point at position `index`.
  */
  std::size_t x(std::size_t index) const
  {
    if (m_lower_scale == NULL)
      return m_x[index];

    // else
    return m_lower_scale->x(index) / 2;
  }
  /*!
    \brief Returns the `y` grid coordinate of the point at position `index`.
  */
  std::size_t y(std::size_t index) const
  {
    if (m_lower_scale == NULL)
      return m_y[index];

    // else
    return m_lower_scale->y(index) / 2;
  }
};
  

}
  
}


#endif // CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H
