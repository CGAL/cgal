// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_SITE_ACCESSORS_H
#define CGAL_VORONOI_DIAGRAM_2_SITE_ACCESSORS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================

template<class S, class DG, class Use_const_ref> struct Site_accessor;
template<class S, class DG, class Use_const_ref> struct Point_accessor;

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

template<class S, class DG, class Use_const_ref>
struct Site_accessor
{
  typedef typename Const_ref_chooser<S,Use_const_ref>::Type  result_type;
  typedef typename DG::Vertex_handle                         Vertex_handle;

  result_type operator()(const Vertex_handle& v) const {
    return v->site();
  }
};

//=========================================================================

template<class P, class DG, class Use_const_ref>
struct Point_accessor
{
  typedef typename Const_ref_chooser<P,Use_const_ref>::Type  result_type;
  typedef typename DG::Vertex_handle                         Vertex_handle;

  result_type operator()(const Vertex_handle& v) const {
    return v->point();
  }
};

//=========================================================================
//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_SITE_ACCESSORS_H
