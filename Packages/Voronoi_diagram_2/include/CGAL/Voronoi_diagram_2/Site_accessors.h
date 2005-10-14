// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_SITE_ACCESSORS_H
#define CGAL_VORONOI_DIAGRAM_2_SITE_ACCESSORS_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class S, class H, class Use_const_ref> struct Site_accessor;
template<class S, class H, class Use_const_ref> struct Point_accessor;

//=========================================================================

template<class T, class Use_const_ref> struct Const_ref_chooser;

template<class T>
struct Const_ref_chooser<T,Tag_true>
{
  typedef const T&  Type;
};

template<class T>
struct Const_ref_chooser<T,Tag_false>
{
  typedef T         Type;
};

//=========================================================================

template<class S, class Vertex_handle, class Use_const_ref>
struct Site_accessor
{
  typedef typename Const_ref_chooser<S,Use_const_ref>::Type  result_type;
  typedef Vertex_handle                                      value_type;
  typedef Arity_tag<1>                                       Arity;

  result_type operator()(const Vertex_handle& v) const {
    return v->site();
  }
};

//=========================================================================

template<class P, class Vertex_handle, class Use_const_ref>
struct Point_accessor
{
  typedef typename Const_ref_chooser<P,Use_const_ref>::Type  result_type;
  typedef Vertex_handle                                      value_type;
  typedef Arity_tag<1>                                       Arity;

  result_type operator()(const Vertex_handle& v) const {
    return v->point();
  }
};

//=========================================================================
//=========================================================================

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_SITE_ACCESSORS_H
