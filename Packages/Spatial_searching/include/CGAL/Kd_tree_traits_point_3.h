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
// file          : include/CGAL/Kd_tree_traits_point_3.h
// package       : ASPAS
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef CGAL_KD_TREE_TRAITS_POINT_3_H
#define CGAL_KD_TREE_TRAITS_POINT_3_H

namespace CGAL {


  template <class K>

  class Kd_tree_traits_point_3 {

  public:
    
    typedef typename K::Cartesian_const_iterator_3 Cartesian_const_iterator;
    typedef typename K::Construct_cartesian_const_iterator_3 Construct_cartesian_const_iterator;
    typedef typename K::Point_3 Point;
    typedef Point** Point_iterator;
    typedef typename K::FT NT;
 
  };

  
} // namespace CGAL
#endif // KD_TREE_TRAITS_POINT_3_H
