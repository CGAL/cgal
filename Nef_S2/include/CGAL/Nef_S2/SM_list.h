// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SM_LIST_H
#define CGAL_SM_LIST_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/In_place_list.h>
#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_2/Object_handle.h>
#include <CGAL/Nef_S2/Generic_handle_map.h>
#include <CGAL/Nef_2/iterator_tools.h>
#include <list>

namespace CGAL {

/*
template <typename HE>
class move_edge_around_svertex {
public:
  void forward(HE& e) const  { e = (e->sprev()->twin()); }
  void backward(HE& e) const { e = (e->twin()->snext()); }
};

template <typename HE>
struct move_edge_around_sface {
  void forward(HE& e)  const { e = (e->snext()); }
  void backward(HE& e) const { e = (e->sprev()); }
};
*/

template < class SVertex>
class SNC_in_place_list_svertex
    : public SVertex, 
      public In_place_list_base<SNC_in_place_list_svertex<SVertex> > {
public:
    typedef SNC_in_place_list_svertex<SVertex> Self;
    //    typedef typename SVertex::SVertex_handle       SVertex_handle;
    //    typedef typename SVertex::SVertex_const_handle SVertex_const_handle;
    SNC_in_place_list_svertex() {}
    SNC_in_place_list_svertex(const SVertex& v)   // down cast
        : SVertex(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SVertex*)this) = ((const SVertex&)v);
        return *this;
    }
};

template < class SHalfedge>
class SNC_in_place_list_shalfedge
    : public SHalfedge, 
      public In_place_list_base<SNC_in_place_list_shalfedge<SHalfedge> > {
public:
    typedef SNC_in_place_list_shalfedge<SHalfedge> Self;
    //    typedef typename SHalfedge::SHalfedge_handle       SHalfedge_handle;
    //    typedef typename SHalfedge::SHalfedge_const_handle SHalfedge_const_handle;
    SNC_in_place_list_shalfedge() {}
    SNC_in_place_list_shalfedge(const SHalfedge& v)   // down cast
        : SHalfedge(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SHalfedge*)this) = ((const SHalfedge&)v);
        return *this;
    }
};

template < class SFace>
class SNC_in_place_list_sface
    : public SFace, 
      public In_place_list_base<SNC_in_place_list_sface<SFace> > {
public:
    typedef SNC_in_place_list_sface<SFace> Self;
    //    typedef typename SFace::SFace_handle       SFace_handle;
    //    typedef typename SFace::SFace_const_handle SFace_const_handle;
    SNC_in_place_list_sface() {}
    SNC_in_place_list_sface(const SFace& v)   // down cast
        : SFace(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SFace*)this) = ((const SFace&)v);
        return *this;
    }
};

} //namespace CGAL

#endif // CGAL_SM_LIST_H
