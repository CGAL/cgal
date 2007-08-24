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
// $URL$
// $Id$
// 
//
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_SHORTEST_PATH_2_H
#define CGAL_SHORTEST_PATH_2_H

#include <CGAL/Visibility_complex_2.h>
#include <functional>
#include <set>
#include <list>
#include <vector>
#include <math.h>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------
namespace Visibility_complex_2_details {
template < class Vc_ >
class Atom 
{
public:
    typedef typename Vc_::Vertex_handle  Vertex_handle;
    typedef typename Vc_::Edge_handle    Edge_handle;
    typedef typename Vc_::Gt::Distance_NT   Distance_NT;
    typedef typename Vc_::Disk_handle Disk_handle;
private:
    Vertex_handle v_;
    Edge_handle   e_;
public:
    Atom() : v_(0) , e_(0) {}
    Atom(Vertex_handle v) : v_(v) , e_(0) {}
    Atom(Edge_handle e) : v_(0) , e_(e) {}
    bool operator==(const Atom& a) const
	{ return (e_ == a.e_ && v_ == a.v_); }
    bool operator!=(const Atom& a) const
	{ return !(*this == a); }

    Edge_handle   edge()        const { return e_; }
    Vertex_handle vertex()      const { return v_; }

    long index()                const 
	{ return (e_ == 0) ? long(v_)          : long(e_);       }
    Distance_NT distance()         const 
	{ return (e_ == 0) ? v_->distance()    : e_->distance(); }
    Distance_NT weight()           const 
	{ return (e_ == 0) ? v_->weight()      : e_->weight();   }
    Atom   prev()               const 
	{ return (e_ == 0) ? v_->prev()        : e_->prev() ;    }
    Vertex_handle next_vertex(Disk_handle s) const 
	{ return (e_ == 0) ? v_->next_vertex(s) : e_->next_vertex(s); }
    Edge_handle next_edge(Disk_handle s)     const 
	{ return (e_ == 0) ? v_->next_edge(s)   : e_->next_edge(s); }
};


template < class Vc_ >
class Sh_edge
    : public Edge_base<Vc_>
{
public:
    typedef typename Vc_::Gt                           Gt;
    typedef typename Vc_::Disk_handle               Disk_handle;
    typedef typename Gt::Distance_NT                      Distance_NT;
    typedef typename Gt::Arc_2                  Arc;
    typedef typename Gt::Point_2                       Point_2;
    typedef Edge_base<Vc_>          Base;
    typedef typename Base::Vertex_handle                        Vertex_handle;
    typedef typename Base::Edge_handle                          Edge_handle;
    typedef Atom<Vc_>                               Atom;

  using Base::inf;
  using Base::sup;
  using Base::object;
  using Base::sign;

  
private:
    Distance_NT      distance_;
    Atom          prev_;
    typename Gt::Length length;

public:
    // -------------------------------------------------------------------------
    Sh_edge() : Base() , distance_(Distance_NT(-1)) { }
    Sh_edge(bool s,Disk_handle p) : Base(s,p), distance_(Distance_NT(-1)) { }
    Sh_edge(Vertex_handle v0 , Vertex_handle v1 , Disk_handle p)
	: Base(v0,v1,p) , distance_(Distance_NT(-1)) { }
    // -------------------------------------------------------------------------
    Distance_NT weight() const { return length(*this,*inf(),*sup()); }
    // -------------------------------------------------------------------------
    Distance_NT distance() const { return distance_; }
    void   set_distance(Distance_NT d) { distance_ = d; }
    // -------------------------------------------------------------------------
    Atom  prev() const { return prev_; }
    void  set_prev(Atom v) { prev_ = v; }
    // -------------------------------------------------------------------------
    Vertex_handle next_vertex(Disk_handle s)
    {
	if (object() != s) {
	    if ((sign()  && sup()->target_object() == object()) ||
	        (!sign() && inf()->target_object() == object())) return 0;
	}
	return (sign()) ? sup() : inf();
    }
    Edge_handle   next_edge(Disk_handle s)
    {
	if (object() != s) {
	    if (( sign() && sup()->is_constraint()) || 
		(!sign() && inf()->is_constraint())) return 0;
	}
	if (sign()) return (sup()->target_object() == object()) ?
			sup()->ccw_target_edge()  : sup()->ccw_source_edge();
	return (inf()->target_object() == object()) ?
	    inf()->cw_target_edge() : inf()->cw_source_edge();
    }
    // -------------------------------------------------------------------------
};



template < class Vc_ >
class Sh_vertex
    : public Vertex_base<Vc_>
{
public:
    typedef typename Vc_::Gt                       Gt;
    typedef typename Vc_::Gt::Distance_NT             Distance_NT;
    typedef typename Vc_::Vertex_handle            Vertex_handle;
    typedef typename Vc_::Edge_handle              Edge_handle;
    typedef typename Vc_::Bitangent_2              Bitangent_2;
    typedef Vertex_base<Vc_>    Base;
    typedef typename Base::Disk_handle                   Disk_handle;
    typedef typename Base::Type                             Type;
    typedef Atom<Vc_>                           Atom;

    using Base::is_constraint;
    using Base::source_object;
    using Base::target_object;
    using Base::is_left_right;
    using Base::is_right_left;
    using Base::is_xx_left;
    using Base::cw_source_edge;
    using Base::cw_target_edge;
    using Base::ccw_source_edge;
    using Base::ccw_target_edge;
private:
    Distance_NT distance_;
    Atom     prev_;
    typename Gt::Length length;

public:

    Sh_vertex() : Base() , distance_(Distance_NT(-1)) { }
    Sh_vertex(Type t , Disk_handle start , Disk_handle finish)
	: Base (t,start,finish) , distance_(Distance_NT(-1)) { }
    Sh_vertex(Edge_handle start , Edge_handle finish)
	: Base(start,finish)    , distance_(Distance_NT(-1)) { }
    Sh_vertex(const Bitangent_2& b) : Base(b) , distance_(Distance_NT(-1)) { }
    Sh_vertex(const Sh_vertex&sibling,bool reverse,Type t)
      : Base(sibling,reverse,t),distance_(Distance_NT(-1)) {}

    Distance_NT weight()   const { return length(*this); }

    Distance_NT distance() const { return distance_; }
    void   set_distance(Distance_NT d) { distance_ = d; }

    Atom prev() const { return prev_; }
    void set_prev(Atom e) { prev_ = e; }

    Vertex_handle next_vertex(Disk_handle /*s*/)  { return 0; }
    Edge_handle   next_edge(Disk_handle s)   {
      if (is_constraint() && s != source_object() && s != target_object()) {
        if (is_left_right() && prev().edge() ==  cw_source_edge()) 
          return 0;
        if (is_right_left() && prev().edge() == ccw_source_edge()) 
          return 0;
      }
      return (is_xx_left()) ? ccw_target_edge() : cw_target_edge(); 
    }
};


class Sh_items : public Items {
public:
    template <class Vc_ >
    struct Edge_wrapper {
	typedef Sh_edge<Vc_>   Edge;
    };
    template <class Vc_>
    struct Vertex_wrapper {
	typedef Sh_vertex<Vc_> Vertex;
    };
};


template <class At>
struct Less_atom {
    bool operator() (const At& a, const At& b) const {
	typedef typename At::Distance_NT Distance_NT;
	if (a.distance() == b.distance()) return (a.index() < b.index());
	if (a.distance() == Distance_NT(-1)) return false;
	if (b.distance() == Distance_NT(-1)) return true;
	return (a.distance() < b.distance());
    }
};

}

typedef Visibility_complex_2_details::Sh_items Shortest_path_2_items;

template <class Gtr_, class Items_=Shortest_path_2_items,
          class FlipTraits=Visibility_complex_2_details::Flip_traits>
class Shortest_path_2 {
  typedef Shortest_path_2<Gtr_,Items_,FlipTraits> Self;
public:
  typedef Gtr_ Gt;
  typedef Items_ Items;
  typedef Visibility_complex_2<Gt,Items,FlipTraits> Visibility_complex_2;
  typedef typename Visibility_complex_2::Vertex Vertex;
  typedef typename Visibility_complex_2::Edge Edge;
  typedef typename Visibility_complex_2::Face Face;

  typedef typename Gt::Distance_NT                     Distance_NT;
  typedef typename Gt::Point_2                      Point_2;
  typedef typename Gt::Bitangent_2                  Bitangent_2;
  typedef typename Gt::Disk         Disk;
  typedef typename Visibility_complex_2::Disk_handle  Disk_handle;
  typedef typename Visibility_complex_2::Edge_handle  Edge_handle;
  typedef typename Visibility_complex_2::Vertex_handle Vertex_handle;

private:
  Visibility_complex_2 vc;
  Edge_handle p,n;
  typedef typename Edge::Atom Atom;
public:

  Shortest_path_2(Visibility_complex_2 V) : vc(V),p(0),n(0) {}
  Shortest_path_2() :p(0),n(0) {};

  void compute_shortest_paths(const Disk& s)
  {
    for (typename Visibility_complex_2::Vertex_iterator i=vc.vertices_begin();
         i!=vc.vertices_end();++i) i->set_distance(-1);
    for (typename Visibility_complex_2::Edge_iterator i=vc.edges_begin();
         i!=vc.edges_end();++i) {
      Edge_handle e=&*i;
      e->set_distance(-1); 
    }
    p=vc.positive_edge(s);
    n=vc.negative_edge(s);
    compute(p);
    compute(n);
  }

  template <class OutputIterator>
  std::pair<Distance_NT,OutputIterator>
  get_path_vertices(const Disk& t, OutputIterator result) {
    if (p&&n) {
      std::list<Vertex_handle> pathp;
      Distance_NT minp_ft = recover_path(p->object(),&t,
                                         std::back_inserter(pathp));
      std::list<Vertex_handle> pathn;
      Distance_NT minn_ft = recover_path(n->object(),&t,
                                         std::back_inserter(pathn));
      // Compare the two paths and return the smallest
      if (minn_ft>=0&&minn_ft<minp_ft) {
        ;
        return std::pair<Distance_NT,OutputIterator>(minn_ft,
                 std::copy(pathn.rbegin(),pathn.rend(),result));
      }
      else if (minp_ft>=0) {
        return std::pair<Distance_NT,OutputIterator>(minp_ft,
                 std::copy(pathp.rbegin(),pathp.rend(),result));
      }
    }
    return std::pair<Distance_NT,OutputIterator>(-1,result);
  }
  template <class OutputIterator>
  std::pair<Distance_NT,OutputIterator>
  get_path_bitangents(const Disk& t, OutputIterator result) {
    std::vector<Vertex_handle> path;
    Distance_NT dist=get_path_vertices(t,std::back_inserter(path)).first;
    for (typename std::vector<Vertex_handle>::const_iterator i=path.begin();
         i!=path.end();++ i) {
      *result=static_cast<const Bitangent_2&>(**i);
      result++;
    }
    return std::pair<Distance_NT,OutputIterator>(dist,result);
  }

private:

  
  // Compute the shortest paths from s
  void compute(Edge_handle s)
  {
    s->set_distance(0);
    s->set_prev(Atom());
    Atom start(s);
    // The priority queue, we push start
    typedef std::set<Atom, Visibility_complex_2_details::Less_atom<Atom> > 
      Queue;
    Queue X;

    X.insert(start);
    // Dijkstra algorithm ------------------------------------------------------
    while (!X.empty()) 
      {
	Atom a = *X.begin(); 
	X.erase(X.begin());

	Distance_NT da = a.distance(); 
	Vertex_handle v = a.next_vertex(s->object());
	Edge_handle   e = a.next_edge(s->object());

	if (v != 0) {
          Distance_NT d = da + v->weight();
          if (v->distance()<0 || d < v->distance()) {
            if (v->distance()>=0) {
              typename Queue::iterator vit = X.find(Atom(v));
              if (vit != X.end()) X.erase(vit);
            }
            v->set_distance(d); v->set_prev(a);
            X.insert(Atom(v));
          }
	}
	if (e != 0) {
          Distance_NT d=da+((e->object()==s->object()&&da==0)?0:e->weight());
          if (e->distance()<0 || d < e->distance()) {
            if (e->distance()>=0) {
              typename Queue::iterator eit = X.find(Atom(e));
              if (eit != X.end()) X.erase(eit);
            }
            e->set_distance(d); e->set_prev(a);
            X.insert(Atom(e));
          }
	}
      }
  }

  template < class OutputIterator>
  Distance_NT
  recover_path(Disk_handle s,Disk_handle t, OutputIterator result)
  {
    // Find edges two edges with opposite sign on t
    Edge_handle ept,emt;
    ept=vc.positive_edge(*t);
    emt=vc.negative_edge(*t);

    // closest edge to s on t
    Edge_handle min = ept;
    Edge_handle f   = ept; 
    do {
      if (min->distance()<0||
          (f->distance() < min->distance()&&f->distance()>=0)) min = f;
      f = f->sup()->ccw_edge(f->object());
    } while (f==ept?(f=emt,true):f!=emt);
    // Shortest path from s to object t
    if (min->distance()<0) return -1;
    recover_path(s,min,result);
    return min->distance();
  }




  template <class OutputIterator>
  OutputIterator
  recover_path(Disk_handle s,Edge_handle t, OutputIterator result)
  {
    Atom finish(t);
    std::list<Vertex_handle> path;
    while (true) {
      if (finish.vertex()!=0) 
        path.push_back(finish.vertex());
      else if (finish.edge()!=0) {
        if (finish.edge()->object()==s)
          return std::copy(path.begin(),path.end(),result);        
      } else
        return result;
      finish=finish.prev();
    }
  }
};

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
