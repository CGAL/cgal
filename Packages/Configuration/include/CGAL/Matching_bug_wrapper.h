// ======================================================================
//
// Copyright (c) 1997-2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Matching_bug_wrapper.h
// package       : Configuration
// maintainer    : Geert-Jan Giezeman <geert@cs.uu.nl>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : Utrecht University
//
// ======================================================================

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
