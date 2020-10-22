// Copyright (c) 2020  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_ORTHTREE_TRAITS_D_H
#define CGAL_ORTHTREE_TRAITS_D_H

namespace CGAL
{

template <typename GeomTraits, typename DimensionTag>
struct Orthtree_traits_d
{
public:

  typedef DimensionTag Dimension;

  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_d Point_d;
  typedef typename GeomTraits::Vector_d Vector_d;
  typedef typename GeomTraits::Sphere_d Sphere_d;
  typedef typename GeomTraits::Cartesian_const_iterator_d Cartesian_const_iterator_d;
  typedef std::array<FT, Dimension::value> Array;

  class Bbox_d
  {
    Point_d m_min, m_max;
  public:

    Bbox_d (const Point_d& min, const Point_d& max)
      : m_min (min), m_max (max)
    { }

    const Point_d& min() { return m_min; }
    const Point_d& max() { return m_max; }
  };

  struct Construct_point_d_from_array
  {
    Point_d operator() (const Array& array) const
    {
      return Point_d ( /* todo */ );
    }
  };
  Construct_point_d_from_array construct_point_d_from_array_object() const
  { return Construct_point_d_from_array(); }

  struct Construct_bbox_d
  {
    Bbox_d operator() (const Array& min,
                       const Array& max) const
    {
      return Bbox_d ( /* todo */ );
    }
  };
  Construct_bbox_d construct_bbox_d_object() const
  { return Construct_bbox_d(); }

};

}

#endif // CGAL_ORTHTREE_TRAITS_D_H
