// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/HDS_items_MSC.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// source        : nef_2d/PM_decorator.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Plane map decorator classes for Mickey Mouse
// ============================================================================

#ifndef CGAL_HDS_ITEMS_MSC_H
#define CGAL_HDS_ITEMS_MSC_H

#include <CGAL/basic.h>
#include <list>

CGAL_BEGIN_NAMESPACE

namespace CGALMSC {

#ifndef CGAL_USE_LEDA
#define LEDA_MEMORY(t) 
#endif

template <typename Refs, typename Traits> class Vertex_;
template <typename Refs, typename Traits> class Halfedge_;
template <typename Refs, typename Traits> class Face_;

template <typename Refs, typename Traits>
class Vertex_  { public:
      typedef Vertex_<Refs,Traits> Base;
      typedef CGAL::Tag_true Supports_vertex_halfedge;
      typedef CGAL::Tag_true Supports_vertex_point;
      typedef typename Refs::Vertex_handle         Vertex_handle;
      typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
      typedef typename Refs::Halfedge_handle       Halfedge_handle;
      typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
      typedef typename Refs::Face_handle           Face_handle;
      typedef typename Refs::Face_const_handle     Face_const_handle;
      typedef typename Refs::Halfedge              Halfedge;
      typedef typename Refs::Face                  Face;
      typedef void*                                GenPtr;

      typedef typename Traits::Point Point;
      typedef typename Traits::Mark  Mark;
      typedef typename std::list<Vertex_handle>::iterator iv_iterator;
    private:
      Halfedge_handle _h;
      Face_handle     _f;
      Point           _p;
      iv_iterator     _ivit;
      Mark            _m;
      GenPtr          _i;
    public:

Vertex_() : _h(),_f(),_ivit(nil_),_m(),_i((GenPtr)0xABCD) {}
Vertex_(const Point& p) : 
  _h(),_f(),_p(p),_ivit(nil_),_m(),_i((GenPtr)0xABCD) {}
     
bool is_isolated() const { return _h == Halfedge_handle(); }
   
Halfedge_handle halfedge() { return _h; }
Halfedge_const_handle halfedge() const { return _h; }
void set_halfedge(Halfedge_handle h) { _h=h; }

Face_handle face() { return _f; }
Face_const_handle face() const { return _f; }
void set_face(Face_handle f) { _f=f; }

Point&  point() { return _p; }
const Point& point() const { return _p; }

Mark& mark() { return _m; }
const Mark& mark() const { return _m; }

GenPtr&       info()       { return _i; }
const GenPtr& info() const { return _i; }

iv_iterator ivit() const { return _ivit; }
void set_ivit(iv_iterator it) { _ivit = it; }

static iv_iterator nil_;

LEDA_MEMORY(Vertex_)

}; // Vertex



template <typename Refs, typename Traits>
class Halfedge_ { public:
      typedef Halfedge_<Refs,Traits> Base;
      typedef Halfedge_<Refs,Traits> Base_base;
      typedef CGAL::Tag_true      Supports_halfedge_prev;
      typedef CGAL::Tag_true      Supports_halfedge_vertex;
      typedef CGAL::Tag_true      Supports_halfedge_face;
      typedef typename Refs::Vertex_handle         Vertex_handle;
      typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
      typedef typename Refs::Halfedge_handle       Halfedge_handle;
      typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
      typedef typename Refs::Face_handle           Face_handle;
      typedef typename Refs::Face_const_handle     Face_const_handle;
      typedef typename Refs::Vertex                Vertex;
      typedef typename Refs::Face                  Face;
      typedef void*                                GenPtr;

      typedef typename std::list<Halfedge_handle>::iterator fc_iterator;
      typedef typename Traits::Mark  Mark;

    protected:

      Halfedge_handle  opp, prv, nxt;
      Vertex_handle    _v;
      Face_handle      _f;    
      fc_iterator      _fcit;
      Mark             _m;
      GenPtr           _i;
    public:
       
Halfedge_() : 
  opp(),prv(),nxt(),_v(),_f(),_fcit(nil_),_m(),_i((GenPtr)0xABCD) {}

Halfedge_handle       opposite()                        { return opp; }
Halfedge_const_handle opposite() const                  { return opp; }

void                  set_opposite(Halfedge_handle h)   { opp = h; }

Halfedge_handle       prev()                            { return prv; }
Halfedge_const_handle prev() const                      { return prv; }
void                  set_prev(Halfedge_handle h)       { prv = h; }

Halfedge_handle       next()                            { return nxt; }
Halfedge_const_handle next() const                      { return nxt; }
void                  set_next(Halfedge_handle h)       { nxt = h; }

Vertex_handle         vertex()                     { return _v; }
Vertex_const_handle   vertex() const               { return _v; }
void                  set_vertex(Vertex_handle v)  { _v = v; }

Face_handle           face()                       { return _f; }
Face_const_handle     face() const                 { return _f; }
void                  set_face(Face_handle f)      { _f = f; }

bool is_border() const { return _f == Face_handle(); }

Mark& mark() { return _m; }
const Mark& mark() const { return _m; }

GenPtr&       info()       { return _i; }
const GenPtr& info() const { return _i; }

fc_iterator fcit() const      { return _fcit; }
void set_fcit(fc_iterator it) { _fcit=it; }

bool is_hole_entry() const 
{ return _fcit != nil_; }

static fc_iterator nil_;
     
LEDA_MEMORY(Halfedge)

};

template <typename Refs, typename Traits>
class Face_ { public:
      typedef Face_<Refs,Traits> Base;
      typedef CGAL::Tag_true                       Supports_face_halfedge;
      typedef typename Refs::Vertex_handle         Vertex_handle;
      typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
      typedef typename Refs::Halfedge_handle       Halfedge_handle;
      typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
      typedef typename Refs::Face_handle           Face_handle;
      typedef typename Refs::Face_const_handle     Face_const_handle;
      typedef typename Refs::Vertex                Vertex;
      typedef typename Refs::Halfedge              Halfedge;
      typedef void*                                GenPtr;

      typedef typename Traits::Mark  Mark;

    class Hole_iterator 
      /*{\Mtypemember iterator for face cycles. Fits the concept 
         |Halfedge_handle|.}*/ 
        : public std::list<Halfedge_handle>::iterator 
    { typedef typename std::list<Halfedge_handle>::iterator Ibase;
      public:
        Hole_iterator() : Ibase() {}
        Hole_iterator(const Ibase& b) : Ibase(b) {}
        Hole_iterator(const Hole_iterator& i) : Ibase(i) {}  
        operator Halfedge_handle() const { return Ibase::operator*(); }
        Halfedge& operator*() { return *(Ibase::operator*()); }
        Halfedge_handle operator->() { return Ibase::operator*(); }
    };

    class Hole_const_iterator : 
      public std::list<Halfedge_handle>::const_iterator 
    { typedef typename std::list<Halfedge_handle>::const_iterator Ibase;
      public:
        Hole_const_iterator() : Ibase() {}
        Hole_const_iterator(const Ibase& b) : Ibase(b) {}
        Hole_const_iterator(const Hole_const_iterator& i) : Ibase(i) {}  
        operator Halfedge_const_handle() const { return Ibase::operator*(); }
        const Halfedge& operator*() { return *(Ibase::operator*()); }
        Halfedge_const_handle operator->() { return Ibase::operator*(); }
    };

    class Isolated_vertex_iterator 
      /*{\Mtypemember iterator for isolated vertices. Fits the concept 
                      |Vertex_handle|.}*/ 
        : public std::list<Vertex_handle>::iterator 
    { typedef typename std::list<Vertex_handle>::iterator Ibase;
      public:
        Isolated_vertex_iterator() : Ibase() {}
        Isolated_vertex_iterator(const Ibase& b) : Ibase(b) {}
        Isolated_vertex_iterator(const Isolated_vertex_iterator& i) 
          : Ibase(i) {}  
        operator Vertex_handle() const { return Ibase::operator*(); }
        Vertex& operator*() { return *(Ibase::operator*()); }
        Vertex_handle operator->() { return Ibase::operator*(); }
    };

    class Isolated_vertex_const_iterator  
      : public std::list<Vertex_handle>::const_iterator 
    { typedef typename std::list<Vertex_handle>::const_iterator Ibase;
      public:
        Isolated_vertex_const_iterator() : Ibase() {}
        Isolated_vertex_const_iterator(const Ibase& b) : Ibase(b) {}
        Isolated_vertex_const_iterator(
          const Isolated_vertex_const_iterator& i) : Ibase(i) {}  
        operator Vertex_const_handle() const { return Ibase::operator*(); }
        const Vertex& operator*() { return *(Ibase::operator*()); }
        Vertex_const_handle operator->() { return Ibase::operator*(); }
    };

    private:
      Halfedge_handle            _e;
      std::list<Halfedge_handle> FC;
      std::list<Vertex_handle>   IV;
      Mark                       _m;
      GenPtr                     _i;
    public:

      Face_() : _e(),_m(),_i((GenPtr)0xABCD) {}
      ~Face_() { FC.clear(); IV.clear(); }

      void store_fc(Halfedge_handle h) 
      { FC.push_back(h); h->set_fcit(--FC.end());  
        CGAL_assertion(h->is_hole_entry()); }

      void remove_fc(Halfedge_handle h) 
      { CGAL_assertion(h->is_hole_entry());
        FC.erase(h->fcit()); h->set_fcit(Halfedge::nil_); }

      void store_iv(Vertex_handle v)
      { IV.push_back(v); v->set_ivit(--IV.end()); }

      void remove_iv(Vertex_handle v)
      { CGAL_assertion(v->is_isolated()); 
        IV.erase(v->ivit()); v->set_ivit(Vertex::nil_); }
      
      Hole_iterator  fc_begin() { return FC.begin(); }
      Hole_iterator  fc_end()   { return FC.end(); }

      Isolated_vertex_iterator  iv_begin() { return IV.begin(); }
      Isolated_vertex_iterator  iv_end()   { return IV.end(); }

      void clear_all_entries()
      { Hole_iterator hit;
        for (hit = fc_begin(); hit!=fc_end(); ++hit) 
          hit->set_fcit(Halfedge::nil_);
        Isolated_vertex_iterator vit;
        for (vit = iv_begin(); vit!=iv_end(); ++vit) 
          vit->set_ivit(Vertex::nil_);
        FC.clear(); IV.clear(); }

      Hole_const_iterator fc_begin() const { return FC.begin(); }
      Hole_const_iterator fc_end() const   { return FC.end(); }
      Isolated_vertex_const_iterator iv_begin() const { return IV.begin(); }
      Isolated_vertex_const_iterator iv_end() const   { return IV.end(); }

      void set_halfedge(Halfedge_handle h)   { _e = h; } 
      Halfedge_handle       halfedge()       { return _e; }
      Halfedge_const_handle halfedge() const { return _e; }

      Mark& mark() { return _m; }
      const Mark& mark() const { return _m; }

      GenPtr&       info()       { return _i; }
      const GenPtr& info() const { return _i; }

  LEDA_MEMORY(Face_)

};


template <typename R,class T> 
typename Vertex_<R,T>::iv_iterator Vertex_<R,T>::nil_;

template <typename R,class T> 
typename Halfedge_<R,T>::fc_iterator Halfedge_<R,T>::nil_;

} // namespace CGALMSC

CGAL_END_NAMESPACE
#endif // CGAL_HDS_ITEMS_MSC_H


