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
// file          : include/CGAL/Search_traits_2.h
// package       : ASPAS
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef CGAL_SEARCH_TRAITS_2_H
#define CGAL_SEARCH_TRAITS_2_H

namespace CGAL {


  template <class K > 

  class Search_traits_2 {

  public:
    typedef typename K::Point_2 Point_d;
    typedef typename K::Iso_rectangle_2 Iso_box_d;
    typedef typename K::Circle_2 Sphere_d;
    typedef typename K::Cartesian_const_iterator_2 Cartesian_const_iterator_d;
    typedef typename K::Construct_cartesian_const_iterator_2 Construct_cartesian_const_iterator_d;

    typedef typename K::Construct_center_2 Construct_center_d;
    typedef typename K::Construct_squared_radius_2 Construct_squared_radius_d;

    typedef typename K::Construct_iso_rectangle_2 Construct_iso_box_d;
    typedef typename K::FT FT;
  };

  
} // namespace CGAL
#endif // CGAL_SEARCH_TRAITS_2_H
