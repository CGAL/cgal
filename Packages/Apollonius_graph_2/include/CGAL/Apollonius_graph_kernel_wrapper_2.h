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
// coordinator   :
//
// ======================================================================



#ifndef CGAL_APOLLONIUS_GRAPH_KERNEL_WRAPPER_2_H
#define CGAL_APOLLONIUS_GRAPH_KERNEL_WRAPPER_2_H

#include <CGAL/Apollonius_graph_short_names_2.h>

#include <CGAL/Apollonius_site_2.h>
#include <CGAL/Cartesian_converter.h>

CGAL_BEGIN_NAMESPACE

template<class Kernel_base_2>
class Apollonius_graph_kernel_wrapper_2 : public Kernel_base_2
{
public:
  typedef CGAL::Apollonius_site_2<Kernel_base_2>  Site_2;
  typedef Kernel_base_2                           Base;

  struct Compare_x_2 : public Kernel_base_2
  {
    typedef Comparison_result   result_type;    
    typedef Arity_tag<2>        Arity;

    Comparison_result operator()(const Site_2& s1,
				 const Site_2& s2) const
    {
      return this->compare_x_2_object()(s1.point(), s2.point());
    }
  };

  struct Compare_y_2 : public Kernel_base_2
  {
    typedef Comparison_result   result_type;
    typedef Arity_tag<2>        Arity;

    Comparison_result operator()(const Site_2& s1,
				 const Site_2& s2) const
    {
      return this->compare_y_2_object()(s1.point(), s2.point());
    }
  };

  struct Orientation_2 : public Kernel_base_2
  {
    typedef Orientation     result_type;
    typedef Arity_tag<3>    Arity;

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


template<class K1, class K2, class Converter >
class Apollonius_graph_cartesian_converter
  : public Converter
{
private:
  typedef typename K2::Site_2                         K2_Site_2;
  typedef typename K2::Point_2                        K2_Point_2;
  typedef Converter                                   Base;
  typedef typename Converter::Number_type_converter   NT_converter;


public:
  bool
  operator()(const bool& b) const {
    return b;
  }

  K2_Point_2
  operator()(const typename K1::Point_2& p) const
  {
    return Base::operator()(p);
  }

  K2_Site_2
  operator()(const typename K1::Site_2& t) const
  {
    NT_converter nt_cv;

    return K2_Site_2( Base::operator()(t.point()),
		      nt_cv(t.weight())
		      );
  }
};


CGAL_END_NAMESPACE


#endif // CGAL_APOLLONIUS_GRAPH_KERNEL_WRAPPER_2_H
