// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// $URL: 
// $Id: 
// 
//
// Author(s)     :     Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF_SNC_INDEXED_ITEMS_H
#define CGAL_NEF_SNC_INDEXED_ITEMS_H
#include <CGAL/Nef_3/Vertex.h>
#include <CGAL/Nef_3/Halfedge.h>
#include <CGAL/Nef_3/Halffacet.h>
#include <CGAL/Nef_3/Volume.h>
#include <CGAL/Nef_3/SHalfedge.h>
#include <CGAL/Nef_3/SHalfloop.h>
#include <CGAL/Nef_3/SFace.h>

CGAL_BEGIN_NAMESPACE

class Index_generator {
  
 public:
  static int get_unique_index()
  {
    static int unique = 0;
    return unique++;
  }
};

class SNC_indexed_items {
 public:
  template <class Refs> class Vertex :    public Vertex_base<Refs> {};
  template <class Refs> class Halffacet : public Halffacet_base<Refs> {};
  template <class Refs> class Volume :    public Volume_base<Refs> {};
  template <class Refs> class SFace :     public SFace_base<Refs> {};

  template <class Refs> class SHalfloop : public SHalfloop_base<Refs> {
    typedef SHalfloop_base<Refs>                    Base;    
    int index;
  public:
    SHalfloop() : Base(), index(0) {}
    SHalfloop(const SHalfloop<Refs>& sl) : Base(sl), index(0) {}
    SHalfloop<Refs>& operator=(const SHalfloop<Refs>& sl) {
      (Base) *this = (Base) sl;
      index = sl.index;
      return *this;
    }

    void set_index(int idx = Index_generator::get_unique_index()) 
    { index = idx; }
    int get_index() const { return index; }
  };

  template <class Refs> class SHalfedge : public SHalfedge_base<Refs> {
    typedef SHalfedge_base<Refs>                    Base;    
    int index;
    int index2;
  public:
    SHalfedge() : Base(), index(0), index2(0) {}
    SHalfedge(const SHalfedge<Refs>& se) 
      : Base(se), index(se.index), index2(se.index2) {}
    SHalfedge<Refs>& operator=(const SHalfedge<Refs>& se) { 
      (Base) *this = (Base) se;
      index = se.index;
      index2 = se.index2;
      return *this;
    }

    void set_index(int idx = Index_generator::get_unique_index()) 
    { index = index2 = idx; }
    int get_index() const { 
      //      CGAL_assertion(index==index2);
      return index; 
    }
    void set_forward_index(int idx)  { index  = idx;}
    void set_backward_index(int idx) { index2 = idx;}
    int get_forward_index() const  { return index;  }
    int get_backward_index() const { return index2; }
    int get_smaller_index() const { return index < index2 ? index : index2; }
  };

  template <class Refs> class SVertex :   public Halfedge_base<Refs> {
    typedef Halfedge_base<Refs>          Base;
    typedef typename Refs::Mark          Mark;
    int index;
  public:
    SVertex() : Base(), index(0) {}
    SVertex(Mark m) : Base(m), index(0) {}
    SVertex(const SVertex<Refs>& sv) : Base(sv) { index = sv.index; }
    SVertex<Refs>& operator=(const SVertex<Refs>& sv) {
      (Base) *this = (Base) sv;
      index = sv.index;
      return *this;
    }

    void set_index(int idx = Index_generator::get_unique_index()) 
    { index = idx; }
    int get_index() const { return index; }
  };
};

CGAL_END_NAMESPACE
#endif // CGAL_NEF_SNC_INDEXED_ITEMS_H
