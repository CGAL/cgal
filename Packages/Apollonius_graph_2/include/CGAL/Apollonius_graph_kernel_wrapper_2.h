// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/Apollonius_graph_kernel_wrapper_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_APOLLONIUS_GRAPH_KERNEL_WRAPPER_2_H
#define CGAL_APOLLONIUS_GRAPH_KERNEL_WRAPPER_2_H


#include <CGAL/Apollonius_site_2.h>
#include <CGAL/Cartesian_converter.h>

CGAL_BEGIN_NAMESPACE

template<class Kernel_base_2>
class Apollonius_graph_kernel_wrapper_2 : public Kernel_base_2
{
public:
  typedef CGAL::Apollonius_site_2<Kernel_base_2>  Site_2;

  struct Compare_x_2 : public Kernel_base_2
  {
    Comparison_result operator()(const Site_2& s1,
				 const Site_2& s2) const
    {
      return this->compare_x_2_object()(s1.point(), s2.point());
    }
  };

  struct Compare_y_2 : public Kernel_base_2
  {
    Comparison_result operator()(const Site_2& s1,
				 const Site_2& s2) const
    {
      return this->compare_y_2_object()(s1.point(), s2.point());
    }
  };

  struct Orientation_2 : public Kernel_base_2
  {
    Orientation operator()(const Site_2& s1,
			   const Site_2& s2,
			   const Site_2& s3) const
    {
      return this->orientation_2_object()(s1.point(),
					  s2.point(),
					  s3.point());
    }
  };

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
  typedef Apollonius_graph_kernel_wrapper_2<_K1> K1;
  typedef Apollonius_graph_kernel_wrapper_2<_K2> K2;

public:
  bool
  operator()(const bool& b) const {
    return b;
  }

  typename K2::Site_2
  operator()(const typename K1::Site_2& wp) const
  {
    Converter c;

    typename K2::Point_2 p(c(wp.x()), c(wp.y()));
    return typename K2::Site_2( p, c(wp.weight()) );
  }
};


CGAL_END_NAMESPACE


#endif // CGAL_APOLLONIUS_GRAPH_KERNEL_WRAPPER_2_H
