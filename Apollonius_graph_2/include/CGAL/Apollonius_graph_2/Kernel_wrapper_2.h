// Copyright (c) 2003,2004,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_2_KERNEL_WRAPPER_2_H
#define CGAL_APOLLONIUS_GRAPH_2_KERNEL_WRAPPER_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Apollonius_graph_2/basic.h>

#include <CGAL/Apollonius_site_2.h>
//#include <CGAL/Cartesian_converter.h>

namespace CGAL {

namespace ApolloniusGraph_2 {

template<class Kernel_base_2>
class Apollonius_graph_kernel_wrapper_2 : public Kernel_base_2
{
public:
  typedef CGAL::Apollonius_site_2<Kernel_base_2>  Site_2;
  typedef Kernel_base_2                           Base;
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
#if defined(_MSC_VER)
  bool operator()(bool b) const {
    return Base::operator()(b);
  }

  K2_Point_2
  operator()(const typename K1::Point_2& p) const
  {
    return Base::operator()(p);
  }
#else
  using Base::operator();
#endif

  K2_Site_2
  operator()(const typename K1::Site_2& t) const
  {
    NT_converter nt_cv;

    return K2_Site_2( Base::operator()(t.point()),
		      nt_cv(t.weight())
		      );
  }
};


} //namespace ApolloniusGraph_2

} //namespace CGAL


#endif // CGAL_APOLLONIUS_GRAPH_2_KERNEL_WRAPPER_2_H
