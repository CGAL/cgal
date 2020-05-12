// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)


#ifndef CGAL_SEARCH_TRAITS_3_H
#define CGAL_SEARCH_TRAITS_3_H

#include <CGAL/license/Spatial_searching.h>

#include <CGAL/Dimension.h>

namespace CGAL {


  template <class K>

  class Search_traits_3 {

  public:

    typedef Dimension_tag<3> Dimension;
    typedef typename K::Cartesian_const_iterator_3 Cartesian_const_iterator_d;
    typedef typename K::Construct_cartesian_const_iterator_3 Construct_cartesian_const_iterator_d;
    typedef typename K::Point_3 Point_d;
    typedef typename K::Iso_cuboid_3 Iso_box_d;
    typedef typename K::Sphere_3 Sphere_d;
    typedef typename K::Construct_iso_cuboid_3 Construct_iso_box_d;

    typedef typename K::Construct_min_vertex_3 Construct_min_vertex_d;
    typedef typename K::Construct_max_vertex_3 Construct_max_vertex_d;
    typedef typename K::Construct_center_3 Construct_center_d;
    typedef typename K::Compute_squared_radius_3 Compute_squared_radius_d;
    typedef typename K::FT FT;

    Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const {
       return Construct_cartesian_const_iterator_d();
    }

  };


} // namespace CGAL
#endif // SEARCH_TRAITS_3_H
