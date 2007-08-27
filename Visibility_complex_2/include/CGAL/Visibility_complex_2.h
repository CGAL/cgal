// Copyright (c) 2001-2004  ENS of Paris (France).
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
// $Source: /CVSROOT/CGAL/Packages/Visibility_complex/include/CGAL/Visibility_complex_2.h,v $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Pierre Angelier, Michel Pocchiola


#ifndef CGAL_VISIBILITY_COMPLEX_2_H
#define CGAL_VISIBILITY_COMPLEX_2_H

#include <CGAL/Visibility_complex_2/Antichain.h>
#include <CGAL/Visibility_complex_2/iterators.h>
#include <CGAL/Visibility_complex_2/Items.h>
#include <CGAL/Visibility_complex_2/Matroidal_items.h>

#include <boost/shared_ptr.hpp>

#include <algorithm>
#include <functional>
#include <vector>

CGAL_BEGIN_NAMESPACE

template < class Tr , 
	   class It   = Visibility_complex_2_items ,
	   class Flip = Visibility_complex_2_details::Flip_traits>
class Visibility_complex_2 
{
public:

  typedef Tr                                           Gt;
  typedef Visibility_complex_2<Tr,It,Flip>             Self;
  typedef Visibility_complex_2_details::Antichain<Tr,It,Flip> Antichain;
  typedef typename Antichain::Disk_handle              Disk_handle;
  typedef typename Antichain::Disk                     Disk;

  typedef typename Antichain::Vertex                   Vertex;
  typedef typename Antichain::Vertex_handle            Vertex_handle;
  typedef typename Antichain::Vertex_const_handle      Vertex_const_handle;
  typedef typename Antichain::Edge                     Edge;
  typedef typename Antichain::Edge_handle              Edge_handle;
  typedef typename Antichain::Edge_const_handle        Edge_const_handle;
  typedef typename Antichain::Face                     Face;
  typedef typename Antichain::Face_handle              Face_handle;
  typedef typename Antichain::Face_const_handle        Face_const_handle;
  typedef typename Antichain::Bitangent_2              Bitangent_2;

  typedef typename Antichain::Ccw_traits               Ccw_traits;
  typedef typename Antichain::Cw_traits                Cw_traits;

private:
  typedef Visibility_complex_2_details::Iterator_base<Self,
            typename std::vector<Vertex_handle>::const_iterator,
            Vertex*> Iterator_base;
  typedef Visibility_complex_2_details::Iterator_base<Self,
            typename std::vector<Vertex_handle>::const_iterator,
            const Vertex*,Iterator_base> Const_iterator_base;
public:
  typedef Visibility_complex_2_details::Vertex_iterator<Iterator_base,
            Vertex,Vertex&,Vertex*>
    Vertex_iterator;
  typedef Visibility_complex_2_details::Vertex_iterator<Const_iterator_base,
            Vertex,const Vertex&,const Vertex*,
            Vertex_iterator>
    Vertex_const_iterator;

  typedef Visibility_complex_2_details::Edge_iterator<Iterator_base,
            Edge,Edge&,Edge*> 
    Edge_iterator;
  typedef Visibility_complex_2_details::Edge_iterator<Const_iterator_base,
            Edge,const Edge&,const Edge*,Edge_iterator>
    Edge_const_iterator;

  typedef Visibility_complex_2_details::Face_iterator<Iterator_base,
            Face,Face&,Face*> 
    Face_iterator;
  typedef Visibility_complex_2_details::Face_iterator<Const_iterator_base,
            Face,const Face&,const Face*,Face_iterator>
    Face_const_iterator;
  // -------------------------------------------------------------------------
  typedef typename Antichain::Linear_sweep_iterator Linear_sweep_iterator;
  typedef typename Antichain::Linear_sweep_const_iterator 
  Linear_sweep_const_iterator;
  typedef typename Antichain::Sweep_iterator        Sweep_iterator;
  typedef typename Antichain::Sweep_const_iterator  Sweep_const_iterator;

  typedef typename Antichain::Constraint_input Constraint_input;
  typedef typename std::vector<Disk>::iterator Disk_iterator;
  typedef typename Antichain::Constraint_iterator Constraint_iterator;
  typedef typename std::vector<Disk>::const_iterator Disk_const_iterator;
  typedef typename Antichain::Constraint_const_iterator
    Constraint_const_iterator;
private:
  class Rep {

    Edge_iterator edges_begin() { return Edge_iterator(vertices.begin()); }
    Edge_iterator edges_end()   { return Edge_iterator(vertices.end());   }
    Face_iterator faces_begin() { return Face_iterator(vertices.begin()); }
    Face_iterator faces_end()   { return Face_iterator(vertices.end());   }
    Vertex_iterator vertices_begin() {
      return Vertex_iterator(vertices.begin());
    }
    Vertex_iterator vertices_end()   {
      return Vertex_iterator(vertices.end());
    }

  public:
    std::vector<Disk> disks;
    Antichain antichain;
    std::vector<Vertex_handle> vertices;
    typedef std::map<Disk_handle,Edge_handle> EM;
    EM pos;
    EM neg;
    Rep() :
            antichain(), vertices() {};
    template <class DiskIterator , class ConstraintIterator> Rep
    (DiskIterator first, DiskIterator last,
     ConstraintIterator  firstc,ConstraintIterator  lastc)  {
      std::copy(first,last,std::back_inserter(disks));
      antichain.~Antichain();
      // This is awfull, but an antichain cannot be copied, and I cannot
      // initialise it in the base clause, because disks would have to be
      // initialised there too, which breaks Sun CC
      new (&antichain) Antichain(false,disks.begin(),disks.end(),firstc,lastc);
      for (typename Antichain::Edge_iterator e=antichain.edges_begin();
           e!=antichain.edges_end();++e) {
        if (e->object()&&e->sign()) {
          pos.insert(typename EM::value_type(e->object(),&*e));
        } else {
          neg.insert(typename EM::value_type(e->object(),&*e));          
        }
      }
    }
    ~Rep() {
      for (Face_iterator i=faces_begin(); i!=faces_end();++i) {
        delete(i.operator ->());
      }
      for (Edge_iterator i=edges_begin(); i!=edges_end();++i) {
        if (i->object()) delete(i.operator ->());
      }
      delete antichain.infinite_face();
      for (typename std::vector<Vertex_handle>::iterator i=vertices.begin();
           i!=vertices.end();++i) {
        if (!(*i)->is_constraint()) {
          delete (*i)->pi();
          delete *i;
        }
      }
      for (Constraint_iterator i=antichain.constraints_begin();
           i!=antichain.constraints_end();++i) {
        delete i->pi()->source_cusp_edge();
        delete i->pi()->target_cusp_edge();
        delete i->pi();
        delete i->source_cusp_edge();
        delete i->target_cusp_edge();
      }
    }
  };
  boost::shared_ptr<Rep> ptr;
public:
  // -------------------------------------------------------------------------
  Visibility_complex_2() { }
  template <class DiskIterator , class ConstraintIterator>
  Visibility_complex_2(DiskIterator first, DiskIterator last,
                       ConstraintIterator  firstc,ConstraintIterator  lastc);
  template <class DiskIterator>
  Visibility_complex_2(DiskIterator first, DiskIterator last) 
  { 
    std::vector<Constraint_input> L;
    *this = Visibility_complex_2(first,last,L.begin(),L.end());
  }
  Visibility_complex_2(const Visibility_complex_2& vc) : ptr(vc.ptr) {}

  int size()             const      { return 2*ptr->vertices.size(); }

  Face_handle infinite_face() {
    return ptr->antichain.infinite_face();
  }
  Face_const_handle infinite_face() const {
    return ptr->antichain.infinite_face();
  }

  bool is_on_convex_hull(Vertex_handle v) const 
  { return ptr->antichain.is_on_convex_hull(v); }

  Vertex_iterator       vertices_begin()       
  { return Vertex_iterator(ptr->vertices.begin());  }
  Vertex_const_iterator vertices_begin() const 
  { return Vertex_iterator(ptr->vertices.begin());  }
  Vertex_iterator       vertices_end()         
  { return Vertex_iterator(ptr->vertices.end());    }
  Vertex_const_iterator vertices_end()   const 
  { return Vertex_iterator(ptr->vertices.end());    }
  // -------------------------------------------------------------------------
  Edge_iterator edges_begin() { return Edge_iterator(ptr->vertices.begin()); }
  Edge_const_iterator edges_begin() const 
  { return Edge_iterator(ptr->vertices.begin()); }
  Edge_iterator edges_end()   { return Edge_iterator(ptr->vertices.end());   }
  Edge_const_iterator edges_end() const  
  { return Edge_iterator(ptr->vertices.end());   }
  // -------------------------------------------------------------------------
  Face_iterator faces_begin() { return Face_iterator(ptr->vertices.begin()); }
  Face_const_iterator faces_begin() const
  { return Face_iterator(ptr->vertices.begin()); }
  Face_iterator faces_end()   { return Face_iterator(ptr->vertices.end());   }
  Face_const_iterator faces_end() const
  { return Face_iterator(ptr->vertices.end());   }

  Disk_iterator disks_begin() {
    return ptr->disks.begin();
  }
  Disk_iterator disks_end() {
    return ptr->disks.end();
  }
  Constraint_iterator constraints_begin() {
    return ptr->antichain.constraints_begin();
  }
  Constraint_iterator constraints_end() {
    return ptr->antichain.constraints_end();
  }
  Disk_const_iterator disks_begin() const {
    return ptr->disks.begin();
  }
  Disk_const_iterator disks_end() const {
    return ptr->disks.end();
  }
  Constraint_const_iterator constraints_begin() const {
    return ptr->antichain.constraints_begin();
  }
  Constraint_const_iterator constraints_end() const {
    return ptr->antichain.constraints_end();
  }

  Edge_handle positive_edge(const Disk&d) {
    return ptr->pos[&d];
  }
  Edge_handle negative_edge(const Disk& d) {
    return ptr->neg[&d];
  }

  Edge_const_handle positive_edge(const Disk&d) const {
    return ptr->pos[&d];
  }
  Edge_const_handle negative_edge(const Disk& d) const {
    return ptr->neg[&d];
  }

  Antichain * antichain() {
    return &(ptr->antichain);
  }
};



template < class Gtr_ , class It , class Flip >
template < class DiskIterator , class ConstraintIterator >
Visibility_complex_2<Gtr_,It,Flip>::Visibility_complex_2(
    DiskIterator first, DiskIterator last ,
    ConstraintIterator  firstc,ConstraintIterator  lastc) :
  ptr(new Rep(first,last,firstc,lastc))
{
  Sweep_iterator v(&(ptr->antichain)) , vend(&(ptr->antichain),0);
  for ( ; v != vend ; ++v) {
    ptr->vertices.push_back(&(*v));
  }
}

template  <class GeomTraits, 
           class Items=Visibility_complex_2_details::Items,
           class Flip=Visibility_complex_2_details::Flip_traits>
class Compute_free_bitangents_2 {

  typedef Visibility_complex_2_details::Antichain<GeomTraits,Items,Flip>
  Antichain;

  typedef typename Antichain::Constraint_iterator   Constraint_iterator;
  typedef typename Antichain::Edge                  Edge;
  typedef typename Antichain::Edge_handle           Edge_handle;
  typedef typename Antichain::Face                  Face;
  typedef typename Antichain::Face_handle           Face_handle;
  typedef typename Antichain::Vertex                Vertex;
  typedef typename Antichain::Vertex_handle         Vertex_handle;
  typedef typename Antichain::Disk                  Disk;
  typedef typename Antichain::Disk_handle           Disk_handle;
public:
  typedef typename Antichain::Bitangent_2           Bitangent_2;
  typedef typename Antichain::Constraint_input      Constraint_input;
  typedef GeomTraits                                Gt;

  template <class DiskIterator, class OutputIterator>
  OutputIterator operator() 
    (DiskIterator first, DiskIterator last,OutputIterator out) const {
    std::vector<Constraint_input> L;
    return operator()(first,last,L.begin(),L.end(),out);
  }

  template <class DiskIterator, class ConstraintIterator, class OutputIterator>
  OutputIterator operator()
    (DiskIterator first, DiskIterator last,
     ConstraintIterator firstc, ConstraintIterator lastc,
     OutputIterator out) const {
    Antichain a(true,first,last,firstc,lastc);
    typedef typename Antichain::Sweep_iterator        Sweep_iterator;
    Sweep_iterator v(&a) , vend(&a,0);
    while (v!=vend) {
      *out=static_cast<Bitangent_2&>(*v);
      ++out;
      ++v;
    }
    struct explore {
      std::set<Edge_handle> edges;
      std::set<Face_handle> faces;
      std::set<Vertex_handle> vertices;

      void meet(Face_handle f) {
        if (!f||!faces.insert(f).second) return;
        meet(f->inf());
        meet(f->sup());
        meet(f->bottom_edge());
        meet(f->top_edge());
      }
      void meet(Edge_handle e) {
        if (!e||!edges.insert(e).second) return;
        meet(e->inf());
        meet(e->sup());
        meet(e->dl());
        meet(e->dr());
        meet(e->ul());
        meet(e->ur());
      }
      void meet(Vertex_handle v) {
        if (!v||!vertices.insert(v).second) return;
        meet(v->raw_pi());
        meet(v->inf());
        meet(v->sup());
        meet(v->ccw_source_edge());
        meet(v->ccw_target_edge());
        meet(v->cw_target_edge());
        meet(v->cw_source_edge());
        meet(v->target_cusp_edge());
        meet(v->source_cusp_edge());
      }
    };
    explore expl;
    expl.meet(&*a.edges_begin());

    for (typename std::set<Face_handle>::iterator i=expl.faces.begin();
         i!=expl.faces.end();i++) {
      if (*i!=a.infinite_face()) delete *i;      
    }
    for (typename std::set<Edge_handle>::iterator i=expl.edges.begin();
         i!=expl.edges.end();i++) {
      if ((*i)->object()) delete *i;
    }
    delete a.infinite_face();
    for (typename std::set<Vertex_handle>::iterator i=expl.vertices.begin();
         i!=expl.vertices.end();i++) {
      if (!(*i)->is_constraint()) delete *i;
    }
    for (Constraint_iterator i=a.constraints_begin();
         i!=a.constraints_end();++i) {
      delete i->pi()->source_cusp_edge();
      delete i->pi()->target_cusp_edge();
      delete i->source_cusp_edge();
      delete i->target_cusp_edge();
      delete i->pi();
    }
    return out;
  }
};


CGAL_END_NAMESPACE

#endif // CGAL_VISIBILITY_COMPLEX_2_H
