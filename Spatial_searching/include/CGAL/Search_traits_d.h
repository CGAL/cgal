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


#ifndef CGAL_SEARCH_TRAITS_D_H
#define CGAL_SEARCH_TRAITS_D_H

#include <CGAL/license/Spatial_searching.h>

#include <CGAL/Dimension.h>

namespace CGAL {


  template <class K,typename D = Dynamic_dimension_tag>
  class Search_traits_d {

  public:
    typedef D Dimension;
    typedef typename K::Cartesian_const_iterator_d Cartesian_const_iterator_d;
    typedef typename K::Construct_cartesian_const_iterator_d Construct_cartesian_const_iterator_d;
    typedef typename K::Construct_iso_box_d Construct_iso_box_d;
    typedef typename K::Construct_min_vertex_d Construct_min_vertex_d;
    typedef typename K::Construct_max_vertex_d Construct_max_vertex_d;
    typedef typename K::Point_d Point_d;
    typedef typename K::Iso_box_d Iso_box_d;
    typedef typename K::Sphere_d Sphere_d;
    typedef typename K::FT FT;

    Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const {
       return Construct_cartesian_const_iterator_d();
    }
  };


} // namespace CGAL
#endif // CGAL_SEARCH_TRAITS_D_H
