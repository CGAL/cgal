// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
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
// file          : include/CGAL/Search_traits_3.h
// package       : ASPAS
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef CGAL_SEARCH_TRAITS_3_H
#define CGAL_SEARCH_TRAITS_3_H

namespace CGAL {


  template <class K>

  class Search_traits_3 {

  public:
    
    typedef typename K::Cartesian_const_iterator_3 Cartesian_const_iterator_d;
    typedef typename K::Construct_cartesian_const_iterator_3 Construct_cartesian_const_iterator_d;
    typedef typename K::Point_3 Point_d;
    typedef typename K::Iso_cuboid_3 Iso_box_d;
    typedef typename K::Sphere_3 Sphere_d;
    typedef typename K::Construct_iso_cuboid_3 Construct_iso_box_d;

    typedef typename K::Construct_vertex_3 Construct_vertex_d;
    typedef typename K::Construct_center_3 Construct_center_d;
    typedef typename K::Compute_squared_radius_3 Compute_squared_radius_d;
    typedef typename K::FT FT;
 
  };

  
} // namespace CGAL
#endif // SEARCH_TRAITS_3_H
