// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : include/CGAL/Kernel_wrapper_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_KERNEL_WRAPPER_2_H
#define CGAL_KERNEL_WRAPPER_2_H


#include <CGAL/Object.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/Cartesian_converter.h>

CGAL_BEGIN_NAMESPACE

template<class Kernel_base_2>
class Kernel_wrapper_2 : public Kernel_base_2
{
public:
  typedef
  CGAL::Weighted_point<typename Kernel_base_2::Point_2,
		       typename Kernel_base_2::RT>       Weighted_point_2;

  typedef CGAL::Object                                   Object_2;
};


template<class Cartesian_converter>
class Extended_cartesian_converter
{};

template<class _K1, class _K2, class Converter >
class Extended_cartesian_converter < Cartesian_converter<_K1,_K2,
  Converter> >
    : public Cartesian_converter<_K1,_K2,Converter>
{
private:
  typedef Kernel_wrapper_2<_K1> K1;
  typedef Kernel_wrapper_2<_K2> K2;

public:
  bool
  operator()(const bool& b) const {
    return b;
  }

  typename K2::Weighted_point_2
  operator()(const typename K1::Weighted_point_2& wp) const
  {
    Converter c;

    typename K2::Point_2 p(c(wp.x()), c(wp.y()));
    return typename K2::Weighted_point_2( p, c(wp.weight()) );
  }
};


CGAL_END_NAMESPACE


#endif // CGAL_KERNEL_WRAPPER_2_H
