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
// file          : include/CGAL/Kd_tree_traits_point.h
// package       : ASPAS
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef CGAL_KD_TREE_TRAITS_POINT_H
#define CGAL_KD_TREE_TRAITS_POINT_H


namespace CGAL {

  template <class NT_, class Point_, class CartesianCoordinateIterator, class ConstructCartesianCoordinateIterator>
  class Kd_tree_traits_point {
    
  public:
    typedef CartesianCoordinateIterator Cartesian_const_iterator;
    typedef ConstructCartesianCoordinateIterator Construct_cartesian_const_iterator;
    typedef Point_ Point;
    typedef NT_ NT;
    
  };

  
} // namespace CGAL
#endif // KD_TREE_TRAITS_POINT_H
