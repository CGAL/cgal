// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-62 $
// release_date  : $CGAL_Date: 2001/05/11 $
//
// file          : include/CGAL/IO/Arr_polyline_traits_iostream.h
// package       : Arrangement (1.82)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================

//#ifdef CGAL_ARR_POLYLINE_TRAITS_H
#ifndef CGAL_ARR_POLYLINE_TRAITS_IOSTREAM_H   
#define CGAL_ARR_POLYLINE_TRAITS_IOSTREAM_H  

#ifndef CGAL_ARR_POLYLINE_TRAITS_H
#include <CGAL/Arr_polyline_traits.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class R, class Container>
std::ostream & 
operator<<(std::ostream & os, 
	   const typename Arr_polyline_traits<R, Container>::Curve_2 & cv)
{
  typedef typename Arr_polyline_traits<R>::Curve_2 Curve_2;
  typedef typename Curve_2::const_iterator         Points_iterator;
  
  os << cv.size() << std::endl;
  for (Points_iterator points_iter = cv.begin(); 
       points_iter != cv.end(); ++points_iter)
    os << " " << *points_iter;

  return os;
}

template <class R, class Container>
std::ostream & 
operator<<(std::ostream & os,  
	   const typename Arr_polyline_traits<R, Container>::X_curve_2 & cv)
{
  typedef typename Arr_polyline_traits<R>::X_curve_2 X_curve_2;
  typedef typename X_curve_2::const_iterator         Points_iterator;
  
  os << cv.size() << std::endl;
  for (Points_iterator points_iter = cv.begin(); 
       points_iter != cv.end(); ++points_iter)
    os << " " << *points_iter;

  return os;
}

template <class R, class Container>
std::istream& 
operator>>(std::istream& in,  
	   typename Arr_polyline_traits<R, Container>::Curve_2 & cv)
{
  typedef typename Arr_polyline_traits<R>::Curve_2   Curve_2;
  typedef typename Curve_2::value_type               Point_2;
  typedef typename Curve_2::size_type                size_type;

  size_type size;

  in >> size;

  for (size_type i = 0; i < size; ++i)
  {
    Point_2 p;
    in >> p;
    cv.push_back(p);  
  }
  
  return in;
}

template <class R, class Container>
std::istream& 
operator>>(std::istream& in, 
	   typename Arr_polyline_traits<R, Container>::X_curve_2 & cv)
{
  typedef typename Arr_polyline_traits<R>::X_curve_2 X_curve_2;
  typedef typename X_curve_2::value_type             Point_2;
  typedef typename X_curve_2::size_type              size_type;

  size_type size;

  in >> size;

  for (size_type i = 0; i < size; ++i)
  {
    Point_2 p;
    in >> p;
    cv.push_back(p);  
  }
  
  return in;
}

CGAL_END_NAMESPACE

#endif
