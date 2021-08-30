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


#ifndef CGAL_KD_TREE_TRAITS_POINT_H
#define CGAL_KD_TREE_TRAITS_POINT_H

#include <CGAL/license/Spatial_searching.h>

#include <CGAL/Dimension.h>


namespace CGAL {
  template <class FT_, class Point, class CartesianCoordinateIterator, class ConstructCartesianCoordinateIterator, typename D = Dynamic_dimension_tag>
  class Search_traits {

  public:

    typedef D Dimension;

    typedef CartesianCoordinateIterator Cartesian_const_iterator_d;
    typedef ConstructCartesianCoordinateIterator Construct_cartesian_const_iterator_d;
    typedef Point Point_d;
    typedef FT_ FT;

    Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const {
       return Construct_cartesian_const_iterator_d();
    }

  };


} // namespace CGAL
#endif //  KD_TREE_TRAITS_POINT_H
