// Copyright (c) 2001,2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_STATIC_FILTERS_H
#define CGAL_STATIC_FILTERS_H

#include <CGAL/basic.h>
#include <CGAL/Static_filters/tools.h>
#include <CGAL/Static_filters/Orientation_2.h>
#include <CGAL/Static_filters/Orientation_3.h>
#include <CGAL/Static_filters/Side_of_oriented_circle_2.h>
#include <CGAL/Static_filters/Side_of_oriented_sphere_3.h>
#include <CGAL/Static_filters/Coplanar_orientation_3.h>
#include <CGAL/Static_filters/Coplanar_side_of_bounded_circle_3.h>

// This traits class gathers optimized predicates written by hand, using
// a few steps of filtering.  It should work if the initial traits has
// cartesian coordinates which fit exactly in doubles.
//
// Purely static filters code has been removed, since it requires additional
// logic and is not plug'n play (requires users providing bounds).
// If it should be provided again, it should probably be separate.

CGAL_BEGIN_NAMESPACE

template < class K_base >
class Static_filters : public K_base
{
public :

  typedef typename K_base::Point_2 Point_2;
  typedef typename K_base::Point_3 Point_3;

  typedef SF_Orientation_2<Point_2>                 Orientation_2;
  typedef SF_Orientation_3<Point_3>                 Orientation_3;
  typedef SF_Side_of_oriented_circle_2<Point_2>     Side_of_oriented_circle_2;
  typedef SF_Side_of_oriented_sphere_3<Point_3>     Side_of_oriented_sphere_3;
  typedef SF_Coplanar_orientation_3<Point_3, Orientation_2>
                                                    Coplanar_orientation_3;
  typedef SF_Side_of_bounded_circle_3<Point_3>
                                            Coplanar_side_of_bounded_circle_3;

  Orientation_2
  orientation_2_object() const
  { return Orientation_2(); }

  Orientation_3
  orientation_3_object() const
  { return Orientation_3(); }

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  { return Side_of_oriented_circle_2(); }

  Side_of_oriented_sphere_3
  side_of_oriented_sphere_3_object() const
  { return Side_of_oriented_sphere_3(); }

  Coplanar_orientation_3
  coplanar_orientation_3_object() const
  { return Coplanar_orientation_3(); }

  Coplanar_side_of_bounded_circle_3
  coplanar_side_of_bounded_circle_3_object() const
  { return Coplanar_side_of_bounded_circle_3(); }
};

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_H
