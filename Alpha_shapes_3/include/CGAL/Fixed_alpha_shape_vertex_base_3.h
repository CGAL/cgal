// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     :  Sebastien Loriot
//

#ifndef CGAL_FIXED_ALPHA_SHAPE_VERTEX_BASE_3_H
#define CGAL_FIXED_ALPHA_SHAPE_VERTEX_BASE_3_H

#include <CGAL/license/Alpha_shapes_3.h>


#include <utility>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/internal/Classification_type.h>

namespace CGAL {

template <class Gt, class Vb = Triangulation_vertex_base_3<Gt> >
class Fixed_alpha_shape_vertex_base_3
  : public Vb
{
public:

  typedef typename Vb::Cell_handle    Cell_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other   Vb2;
    typedef Fixed_alpha_shape_vertex_base_3<Gt, Vb2>              Other;
  };

  typedef typename Vb::Point Point;
  
private:
  typedef internal::Classification_type Classification_type;  
  Classification_type status_;

  bool is_on_chull_;

public:

  Fixed_alpha_shape_vertex_base_3()    
    : Vb(),is_on_chull_(false) {}
  
  Fixed_alpha_shape_vertex_base_3(const Point& p)
    : Vb(p),is_on_chull_(false) {}
  
  Fixed_alpha_shape_vertex_base_3(const Point& p, Cell_handle c)
    : Vb(p, c),is_on_chull_(false) {}

  Classification_type get_classification_type() { return status_;}
  void set_classification_type(Classification_type status) {status_=status;}
  
  void is_on_chull(bool b){is_on_chull_=b;};
  bool is_on_chull(){return is_on_chull_;}
    
};

} //namespace CGAL

#endif // CGAL_FIXED_ALPHA_SHAPE_VERTEX_BASE_3_H
