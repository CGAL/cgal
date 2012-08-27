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

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 83
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

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
  template <class Refs> class Volume :    public Volume_base<Refs> {};
  template <class Refs> class SFace :     public SFace_base<Refs> {};

  template <class Refs> class Halffacet : public Halffacet_base<Refs> {
    typedef Halffacet_base<Refs>    Base;
    typedef typename Refs::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;
    typedef typename Refs::SHalfedge_const_handle
      SHalfedge_const_handle;
    typedef typename Refs::SHalfloop_const_handle
      SHalfloop_const_handle;
    typedef typename Refs::Halffacet_cycle_const_iterator
      Halffacet_cycle_const_iterator;
    typedef typename Refs::Plane_3  Plane_3;
    typedef typename Refs::Mark     Mark;

  public:
    Halffacet() : Base() {}
    Halffacet(const Plane_3& h, Mark m) : Base(h, m) {}

    bool is_valid( bool verb = false, int level = 0) const {
      bool valid = Base::is_valid(verb, level);

      Halffacet_cycle_const_iterator 
	fci(this->facet_cycles_begin());
      if(!fci.is_shalfedge()) return false;
      SHalfedge_const_handle set(fci);
      int index = set->get_index();
      for(; fci != this->facet_cycles_end(); ++fci) {
	if(fci.is_shalfedge()) {
	  SHalfedge_const_handle se(fci);
	  SHalfedge_around_facet_const_circulator 
	    sfc(se), send(sfc);
	  do {	    
	    valid = valid && sfc->get_index() == index;
	    ++sfc;
	  } while(sfc != send);
	} else if(fci.is_shalfloop()) {
	  SHalfloop_const_handle sl(fci);
	  valid = valid && sl->get_index() == index;
	} else
	  return false;
      }
      return valid;
    }
  };

  template <class Refs> class SHalfloop : public SHalfloop_base<Refs> {
    typedef SHalfloop_base<Refs>                    Base;
    typedef typename Refs::Halffacet_const_handle Halffacet_const_handle;
    int index;
    Halffacet_const_handle ifacet;
    bool init_ifacet;
  public:
    SHalfloop() : Base(), index(0), init_ifacet(false) {}
    SHalfloop(const SHalfloop<Refs>& sl) 
      : Base(sl), index(0), 
      ifacet(sl.ifacet), init_ifacet(sl.init_ifacet) {}
    SHalfloop<Refs>& operator=(const SHalfloop<Refs>& sl) {
      (Base&) *this = (Base) sl;
      index = sl.index;
      ifacet = sl.ifacet;
      init_ifacet = sl.init_ifacet;
      return *this;
    }

    void set_index(int idx = Index_generator::get_unique_index()) 
    { index = idx; }
    int get_index() const { return index; }
    Halffacet_const_handle get_index_facet() const { 
      if(init_ifacet)
	return ifacet;
      return this->facet();
    }
    void set_index_facet(Halffacet_const_handle f) { 
      ifacet = f;
      init_ifacet = true;
    }
  };

  template <class Refs> class SHalfedge : public SHalfedge_base<Refs> {
    typedef SHalfedge_base<Refs>                    Base;
    typedef typename Refs::Halffacet_const_handle Halffacet_const_handle;
    int index;
    int index2;
    Halffacet_const_handle ifacet;
    bool init_ifacet;
  public:
    SHalfedge() : Base(), index(0), index2(0), init_ifacet(false) {}
    SHalfedge(const SHalfedge<Refs>& se)
      : Base(se), index(se.index), index2(se.index2), 
      ifacet(se.ifacet), init_ifacet(se.init_ifacet) {}
    SHalfedge<Refs>& operator=(const SHalfedge<Refs>& se) { 
      (Base&) *this = (Base) se;
      index = se.index;
      index2 = se.index2;
      ifacet = se.ifacet;
      init_ifacet = se.init_ifacet;
      return *this;
    }

    void set_index(int idx = Index_generator::get_unique_index()) 
    { index = index2 = idx; }
    int get_index() const { 
      return index; 
    }
    void set_forward_index(int idx)  { index  = idx;}
    void set_backward_index(int idx) { index2 = idx;}
    int get_forward_index() { return index;  }
    int get_backward_index() { return index2; }
    int get_smaller_index() { return index < index2 ? index : index2; }
    Halffacet_const_handle get_index_facet() const { 
      if(init_ifacet) 
	return ifacet;
      return this->facet();
    }
    void set_index_facet(Halffacet_const_handle f) { 
      ifacet = f;
      init_ifacet = true;
    }
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
      (Base&) *this = (Base) sv;
      index = sv.index;
      return *this;
    }

    void set_index(int idx = Index_generator::get_unique_index()) 
    { index = idx; }
    int get_index() const { return index; }
  };
};

} //namespace CGAL
#endif // CGAL_NEF_SNC_INDEXED_ITEMS_H
