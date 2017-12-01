// Copyright (c) 1997, 2012 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Tran Kai Frank DA

#ifndef CGAL_ALPHA_SHAPE_CELL_BASE_3_H
#define CGAL_ALPHA_SHAPE_CELL_BASE_3_H

#include <CGAL/license/Alpha_shapes_3.h>

#include <vector>
#include <CGAL/Compact_container.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/internal/Lazy_alpha_nt_3.h>
#include <CGAL/Default.h>

namespace CGAL {

template < class NT_>
class  Alpha_status 
: public Compact_container_base
{
  bool _is_Gabriel;
  bool _is_on_chull;
  NT_ _alpha_min;
  NT_ _alpha_mid;
  NT_ _alpha_max;

public:
  typedef NT_ NT;
  Alpha_status() : _is_Gabriel(false), _is_on_chull(false) {}
  void set_alpha_min(NT alpha) {_alpha_min = alpha;}
  void set_alpha_mid(NT alpha) {_alpha_mid = alpha;}
  void set_alpha_max(NT alpha) {_alpha_max = alpha;}
  void set_is_Gabriel(bool yesorno) { _is_Gabriel = yesorno;}
  void set_is_on_chull(bool yesorno) {_is_on_chull = yesorno;}
  NT alpha_min() const { return _alpha_min;}
  NT alpha_mid() const { return _alpha_mid;}
  NT alpha_max() const { return _alpha_max;}
  bool is_Gabriel() const {return _is_Gabriel;}
  bool is_on_chull() const {return  _is_on_chull;}
};

template < class Gt,
           class Cb_ = Default,
           class ExactAlphaComparisonTag = Tag_false,
           class Weighted_tag = Tag_false >
class Alpha_shape_cell_base_3
  : public Default::Get<Cb_, Triangulation_cell_base_3<Gt> >::type
{
  typedef typename Default::Get<Cb_, Triangulation_cell_base_3<Gt> >::type Cb;
public:
  typedef typename Cb::Vertex_handle   Vertex_handle;
  typedef typename Cb::Cell_handle     Cell_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other   Cb2;
    typedef Alpha_shape_cell_base_3<Gt, Cb2,ExactAlphaComparisonTag,Weighted_tag>                Other;
  };

  typedef typename internal::Alpha_nt_selector_3<
    Gt,ExactAlphaComparisonTag,Weighted_tag>::Type_of_alpha  NT;
  typedef CGAL::Alpha_status<NT>                             Alpha_status;
  typedef Compact_container<Alpha_status>                    Alpha_status_container;
  typedef typename Alpha_status_container::const_iterator    Alpha_status_const_iterator;
  typedef typename Alpha_status_container::iterator          Alpha_status_iterator;

private:
  Alpha_status_iterator facet_status[4];
  NT A;

public:
  Alpha_shape_cell_base_3() 
    : Cb() {}
  
  Alpha_shape_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {}
  
  Alpha_shape_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2, Vertex_handle v3,
                          Cell_handle n0, Cell_handle n1,
                          Cell_handle n2, Cell_handle n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3) {}

  NT get_alpha() const { return A; }
  void set_alpha(const NT & AA) { A = AA;}

  Alpha_status_iterator get_facet_status(int i) {return facet_status[i]; }
  Alpha_status_const_iterator get_facet_status(int i) const { return facet_status[i]; }

  void set_facet_status(int i, Alpha_status_iterator as) { facet_status[i]= as; }
};

} // namespace CGAL

#endif // CGAL_ALPHA_SHAPE_CELL_BASE_3_H
