// Copyright (c) 2020  GF (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_OCTREE_TRAITS_3_H
#define CGAL_OCTREE_TRAITS_3_H

namespace CGAL
{

template <typename GeomTraits>
struct Octree_traits_2
{
public:

  typedef Dimension_tag<2> Dimension;
  typedef Bbox_2 Bbox_d;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_2 Point_d;
  typedef typename GeomTraits::Vector_2 Vector_d;
  typedef typename GeomTraits::Iso_rectangle_2 Iso_box_d;
  typedef typename GeomTraits::Circle_2 Sphere_d;
  typedef typename GeomTraits::Aff_transformation_2 Aff_transformation_d;
  typedef typename GeomTraits::Cartesian_const_iterator_2 Cartesian_const_iterator_d;

  struct Construct_point_d_from_array
  {
    Point_d operator() (const std::array<FT, Dimension::value>& array) const
    {
      return Point_d (array[0], array[1]);
    }
  };
  Construct_point_d_from_array construct_point_d_from_array_object() const
  { return Construct_point_d_from_array(); }

  struct Construct_bbox_d
  {
    Bbox_d operator() (const std::array<FT, Dimension::value>& min,
                       const std::array<FT, Dimension::value>& max) const
    {
      return Bbox_d (min[0], min[1], max[0], max[1]);
    }
  };
  Construct_bbox_d construct_bbox_d_object() const
  { return Construct_bbox_d(); }

  struct Construct_iso_bounding_box_d
  {
    template <typename InputIterator>
    Iso_box_d operator() (InputIterator begin, InputIterator end) const
    {
      return CGAL::bounding_box (begin, end);
    }
  };
  Construct_iso_bounding_box_d construct_iso_bounding_box_d_object() const
  { return Construct_iso_bounding_box_d(); }

};

}

#endif // CGAL_OCTREE_TRAITS_2_H
