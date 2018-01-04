// Copyright (c) 1997, 2012  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Tran Kai Frank DA

#ifndef CGAL_ALPHA_SHAPE_VERTEX_BASE_3_H
#define CGAL_ALPHA_SHAPE_VERTEX_BASE_3_H

#include <CGAL/license/Alpha_shapes_3.h>


#include <utility>
#include <CGAL/Compact_container.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Default.h>

namespace CGAL {

template <class Gt, class Vb_ = Default, class ExactAlphaComparisonTag=Tag_false,class Weighted_tag=Tag_false  >
class Alpha_shape_vertex_base_3
  : public Default::Get<Vb_, Triangulation_vertex_base_3<Gt> >::type
{
  typedef typename Default::Get<Vb_, Triangulation_vertex_base_3<Gt> >::type Vb;
public:

  typedef typename Vb::Cell_handle    Cell_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other   Vb2;
    typedef Alpha_shape_vertex_base_3<Gt, Vb2,ExactAlphaComparisonTag,Weighted_tag>              Other;
  };

  typedef typename Vb::Point Point;
  typedef typename internal::Alpha_nt_selector_3<Gt,ExactAlphaComparisonTag,Weighted_tag>::Type_of_alpha  NT;
  typedef CGAL::Alpha_status<NT>     Alpha_status;
  typedef Compact_container<Alpha_status>   Alpha_status_container;
  typedef typename Alpha_status_container::const_iterator 
                                            Alpha_status_const_iterator;
  typedef typename Alpha_status_container::iterator 
                                            Alpha_status_iterator;
  
private:
  Alpha_status  _as;


public:

  Alpha_shape_vertex_base_3()    
    : Vb() {}
  
  Alpha_shape_vertex_base_3(const Point& p)
    : Vb(p) {}
  
  Alpha_shape_vertex_base_3(const Point& p, Cell_handle c)
    : Vb(p, c) {}

  Alpha_status*  get_alpha_status() { return &_as;}
  void set_alpha_status(Alpha_status_iterator as) {_as= as;}
};

} //namespace CGAL

#endif // CGAL_ALPHA_SHAPE_VERTEX_BASE_3_H
