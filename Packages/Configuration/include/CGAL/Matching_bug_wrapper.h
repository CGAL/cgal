// Copyright (c) 1997-2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
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
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_MATCHING_BUG_WRAPPER_H
#define CGAL_MATCHING_BUG_WRAPPER_H

CGAL_BEGIN_NAMESPACE

namespace CGALi {

  // This class is used as a wrapper for template arguments which are
  // a model for Kernel: it just replicates the type information.
  // This is used to make compilers having MATCHING_BUG_4 accept 
  // overloaded functions with arguments of type K::some_type.
  
  template < class K >
  struct Matching_bug_wrapper {
    typedef typename K::Point_2           Point_2;
    typedef typename K::Vector_2          Vector_2;
    typedef typename K::Direction_2       Direction_2;
    typedef typename K::Line_2            Line_2;
    typedef typename K::Ray_2             Ray_2;
    typedef typename K::Segment_2         Segment_2;
    typedef typename K::Triangle_2        Triangle_2;
    typedef typename K::Iso_rectangle_2   Iso_rectangle_2; 
    typedef typename K::Circle_2          Circle_2; 
    // this is not yet part of the kernel models:
    //typedef typename K::Weighted_point_2  Weighted_point_2;   
    typedef typename K::Object_2          Object_2; 

    typedef typename K::Point_3           Point_3;
    typedef typename K::Vector_3          Vector_3;
    typedef typename K::Direction_3       Direction_3;
    typedef typename K::Line_3            Line_3;
    typedef typename K::Ray_3             Ray_3;
    typedef typename K::Segment_3         Segment_3;
    typedef typename K::Triangle_3        Triangle_3;
    typedef typename K::Iso_cuboid_3      Iso_cuboid_3; 
    typedef typename K::Sphere_3          Sphere_3; 
    typedef typename K::Plane_3           Plane_3;
    typedef typename K::Tetrahedron_3     Tetrahedron_3;
    // this is not yet part of the kernel models:
    //typedef typename K::Weighted_point_3  Weighted_point_3;   
    typedef typename K::Object_3          Object_3; 
  };

}

CGAL_END_NAMESPACE

#endif // CGAL_MATCHING_BUG_WRAPPER_H
