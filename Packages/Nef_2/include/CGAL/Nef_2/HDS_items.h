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
// file          : include/CGAL/Nef_2/HDS_items.h
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
// implementation: Extended item classes
// ============================================================================

#ifndef CGAL_HDS_ITEMS_H
#define CGAL_HDS_ITEMS_H

#include <CGAL/basic.h>
#include <list>

template <typename Refs >
struct Halfedge__base {
  typedef typename Refs::Halfedge_handle Halfedge_handle;
  typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
protected:
  Halfedge_handle opp; 
public:
  Halfedge_handle       opposite()                        { return opp; }
  Halfedge_const_handle opposite() const                  { return opp; }
  void                  set_opposite(Halfedge_handle h)   { opp = h; }
};


#ifndef CGAL_USE_LEDA
#define LEDA_MEMORY(t) 
#endif

struct HDS_items {
  /*{\Moptions outfile=epm_vertex.man }*/
  /*{\Moptions constref=yes}*/
  /*{\Moptions print_title=yes }*/ 
  /*{\Moptions section=subsection}*/
  /*{\Manpage{Vertex}{}{The HDS vertex base class}{v}}*/    

  template <typename Refs, typename Traits>
  class Vertex_wrapper { public:
    typedef typename Traits::Point Point;
    class Vertex {
    public:
      typedef Refs     HalfedgeDS;
      typedef Vertex   Base;
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

  /*{\Mtypes 3}*/

      typedef typename Traits::Point Point;
      /*{\Mtypemember geometric embedding}*/ 

      typedef typename Traits::Mark  Mark;
      /*{\Mtypemember information}*/ 

      typedef typename std::list<Vertex_handle>::iterator iv_iterator;
    private:
      Halfedge_handle _h;
      Face_handle     _f;
      Point           _p;
      iv_iterator     _ivit;
      Mark            _m;
      GenPtr          _i;
    public:

  /*{\Mcreation 3}*/ 
      Vertex() : 
        _h(),_f(),_ivit(nil_),_m(),_i((GenPtr)0xABCD) {}
      /*{\Mcreate constructs an uninitialized vertex concerning embedding 
      and mark. All links are initialized by their default value.}*/   

      Vertex(const Point& p) : 
        _h(),_f(),_p(p),_ivit(nil_),_m(),_i((GenPtr)0xABCD) {}
      /*{\Mcreate constructs a vertex with embedding |p| and mark |m|.
         All links are initialized by their default value.}*/   

  /*{\Moperations 3 4}*/
     
      bool is_isolated() const 
      /*{\Mop returns true iff |\Mvar| is isolated, else false.}*/
      { return _h == Halfedge_handle(); }
   
      Halfedge_handle halfedge() { return _h; }
      /*{\Mop returns an incident halfedge. \precond |!is_isolated()|.}*/
      Halfedge_const_handle halfedge() const { return _h; }

      void set_halfedge(Halfedge_handle h) { _h=h; }
      /*{\Mop makes |h| the entry point into the adjacency cycle of 
      |\Mvar|.}*/

      Face_handle face() { return _f; }
      /*{\Mop returns the incident face if |is_isolated()|.}*/
      Face_const_handle face() const { return _f; }

      void set_face(Face_handle f) { _f=f; }
      /*{\Mop makes |f| the incident face of |\Mvar|.}*/

      Point& point() { return _p; }
      /*{\Mop returns the embedding point of |\Mvar|.}*/
      const Point& point() const { return _p; }

      Mark& mark() { return _m; }
      /*{\Mop returns the mark of |\Mvar|.}*/
      const Mark& mark() const { return _m; }

      GenPtr& info() { return _i; }
      /*{\Mop returns a generic information slot of |\Mvar|.}*/
      const GenPtr& info() const { return _i; }

      iv_iterator ivit() const { return _ivit; }
      void set_ivit(iv_iterator it) { _ivit = it; }
      static iv_iterator nil_;
      /* stl iterators have default construction but are only equal 
      comparable when copy constructed, what a mess in the specification */

      LEDA_MEMORY(Vertex)
    };
  };


  /*{\Moptions outfile=epm_halfedge.man}*/
  /*{\Moptions constref=yes}*/
  /*{\Moptions print_title=yes }*/ 
  /*{\Moptions section=subsection}*/
  /*{\Manpage{Halfedge}{}{The HDS halfedge base class}{e}}*/ 

  template <typename Refs, typename Traits>
  class Halfedge_wrapper { public:
    struct Halfedge {
    public:
      typedef Refs                HalfedgeDS;
      typedef Halfedge            Base;
      typedef Halfedge            Base_base;
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

  /*{\Mtypes 3}*/
      typedef typename Traits::Mark  Mark;
      /*{\Mtypemember information}*/ 

    protected:

      Halfedge_handle  opp, prv, nxt;
      Vertex_handle    _v;
      Face_handle      _f;    
      fc_iterator      _fcit;
      Mark             _m;
      GenPtr           _i;
    public:
       
  /*{\Mcreation 3}*/ 
      Halfedge() : 
        opp(),prv(),nxt(),_v(),_f(),_fcit(nil_),_m(),_i((GenPtr)0xABCD) {}
      /*{\Mcreate constructs an uninitialized halfedge concerning embedding 
      and mark. All links are initialized by their default value.}*/   

  /*{\Moperations 3 4}*/

      Halfedge_handle       opposite()                        { return opp; }
      /*{\Mop returns the twin of |\Mvar|.}*/
      Halfedge_const_handle opposite() const                  { return opp; }

      void                  set_opposite(Halfedge_handle h)   { opp = h; }
      /*{\Mop makes |h| the twin of |\Mvar|.}*/

      Halfedge_handle       prev()                            { return prv; }
      /*{\Mop returns the previous edge of the face cycle of |\Mvar|.}*/
      Halfedge_const_handle prev() const                      { return prv; }

      void                  set_prev(Halfedge_handle h)       { prv = h; }
      /*{\Mop makes |h| the previous edge in the face cycle of |\Mvar|.}*/

      Halfedge_handle       next()                            { return nxt; }
      /*{\Mop returns the next edge of the face cycle of |\Mvar|.}*/
      Halfedge_const_handle next() const                      { return nxt; }

      void                  set_next(Halfedge_handle h)       { nxt = h; }
      /*{\Mop makes |h| the next edge in the face cycle of |\Mvar|.}*/

      Vertex_handle         vertex()                     { return _v; }
      /*{\Mop returns the vertex incident to the halfedge |\Mvar|.}*/
      Vertex_const_handle   vertex() const               { return _v; }

      void                  set_vertex(Vertex_handle v)  { _v = v; }
      /*{\Mop makes |v| the vertex incident to |\Mvar|.}*/

      Face_handle           face()                       { return _f; }
      /*{\Mop returns the face incident to the halfedge |\Mvar|.}*/
      Face_const_handle     face() const                 { return _f; }

      void                  set_face(Face_handle f)      { _f = f; }
      /*{\Mop makes |f| the face incident to |\Mvar|.}*/

      bool is_border() const { return _f == Face_handle(); }

      Mark& mark() { return _m; }
      /*{\Mop returns the mark of |\Mvar|.}*/
      const Mark& mark() const { return _m; }

      GenPtr&       info()       { return _i; }
      /*{\Mop returns a generic information slot of |\Mvar|.}*/
      const GenPtr& info() const { return _i; }

      fc_iterator fcit() const      { return _fcit; }
      void set_fcit(fc_iterator it) { _fcit=it; }

      bool is_hole_entry() const 
      /*{\Mop returns true iff |\Mvar| is entry point into a hole face
         cycle of |\Mvar.face()|.}*/
      { return _fcit != nil_; }

      static fc_iterator nil_;
      /* stl iterators have default construction but are only equal comparable
         when copy constructed, what a mess in the specification */
      
      LEDA_MEMORY(Halfedge)
    };
  };

  /*{\Moptions outfile=epm_face.man}*/
  /*{\Moptions constref=yes}*/
  /*{\Moptions print_title=yes }*/ 
  /*{\Moptions section=subsection}*/
  /*{\Manpage{Face}{}{The HDS face base class}{f}}*/    

  template <typename Refs, typename Traits>
  class Face_wrapper { public:
    class Face {
    public:
      typedef CGAL::Tag_true                       Supports_face_halfedge;
      typedef Refs  HalfedgeDS;
      typedef Face  Base;

      typedef typename Refs::Vertex_handle         Vertex_handle;
      typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
      typedef typename Refs::Halfedge_handle       Halfedge_handle;
      typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
      typedef typename Refs::Face_handle           Face_handle;
      typedef typename Refs::Face_const_handle     Face_const_handle;
      typedef typename Refs::Vertex                Vertex;
      typedef typename Refs::Halfedge              Halfedge;
      typedef void*                                GenPtr;

  /*{\Mtypes 3}*/

      typedef typename Traits::Mark  Mark;
      /*{\Mtypemember mark information}*/ 

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

    /*{\Mtext |Hole_const_iterator| and |Isolated_vertex_const_iterator| 
              are the non mutable versions.}*/

    private:
      Halfedge_handle            _e;
      std::list<Halfedge_handle> FC;
      std::list<Vertex_handle>   IV;
      Mark                       _m;
      GenPtr                     _i;
    public:

  /*{\Mcreation 4}*/ 
      Face() : _e(),_m(),_i((GenPtr)0xABCD) {}
      /*{\Mcreate constructs an uninitialized face with undefined mark,
        empty face cycle list, and empty isolated vertices list.}*/   

      ~Face() { FC.clear(); IV.clear(); }

  /*{\Moperations 2.5 3}*/

      void store_fc(Halfedge_handle h) 
      /*{\Mop stores halfedge |h| as an entry into a face cycle of |\Mvar|.
          Postcondition: |h->is_hole_entry()|.}*/
      { FC.push_back(h); h->set_fcit(--FC.end());  
        CGAL_assertion(h->is_hole_entry()); }

      void remove_fc(Halfedge_handle h) 
      /*{\Mop removes halfedge |h| as an entry into a face cycle of |\Mvar|.
          \precond |h->is_hole_entry()| and |h|  is stored in the
          face cycle list of |\Mvar|.
          Postcondition: |!h->is_hole_entry()|.}*/
      { CGAL_assertion(h->is_hole_entry());
        FC.erase(h->fcit()); h->set_fcit(Halfedge::nil_); }

      void store_iv(Vertex_handle v)
      /*{\Mop stores vertex |v| as an isolated vertex of |\Mvar|.}*/
      { IV.push_back(v); v->set_ivit(--IV.end()); }

      void remove_iv(Vertex_handle v)
      /*{\Mop removes vertex |v| as an isolated vertex of |\Mvar|.
          \precond |v->is_isolated()| and |v| is stored in the
          isolated vertices list of |\Mvar|. 
          Postcondition: |!v->is_isolated()|.}*/
      { CGAL_assertion(v->is_isolated()); 
        IV.erase(v->ivit()); v->set_ivit(Vertex::nil_); }
      
      /*{\Mtext\setopdims{4cm}{0cm}}*/

      Hole_iterator  fc_begin() { return FC.begin(); }
      /*{\Mop}*/

      Hole_iterator  fc_end()   { return FC.end(); }
      /*{\Mop the iterator range |[fc_begin(),fc_end())| spans the set of 
          interior face cycles.}*/

      Isolated_vertex_iterator  iv_begin() { return IV.begin(); }
      /*{\Mop}*/

      Isolated_vertex_iterator  iv_end()   { return IV.end(); }
      /*{\Mop the iterator range |[iv_begin(),iv_end())| spans the set of 
        isolated vertices.}*/

      void clear_all_entries()
      { Hole_iterator hit;
        for (hit = fc_begin(); hit!=fc_end(); ++hit) 
          hit->set_fcit(Halfedge::nil_);
        Isolated_vertex_iterator vit;
        for (vit = iv_begin(); vit!=iv_end(); ++vit) 
          vit->set_ivit(Vertex::nil_);
        FC.clear(); IV.clear(); }

  /*{\Mtext There are the same iterator ranges defined for the const
  iterators |Hole_const_iterator|, |Isolated_vertex_const_iterator|.
  \restoreopdims}*/

      Hole_const_iterator fc_begin() const { return FC.begin(); }
      Hole_const_iterator fc_end() const   { return FC.end(); }
      Isolated_vertex_const_iterator iv_begin() const { return IV.begin(); }
      Isolated_vertex_const_iterator iv_end() const   { return IV.end(); }

      void set_halfedge(Halfedge_handle h)   { _e = h; } 
      /*{\Mop makes |h| the entry edge into the outer face cycle.}*/
      Halfedge_handle       halfedge()       { return _e; }
      /*{\Mop returns a halfedge in the outer face cycle.}*/
      Halfedge_const_handle halfedge() const { return _e; }

      Mark& mark() { return _m; }
      /*{\Mop returns the mark of |\Mvar|.}*/
      const Mark& mark() const { return _m; }

      GenPtr&       info()       { return _i; }
      /*{\Mop returns a generic information slot of |\Mvar|.}*/
      const GenPtr& info() const { return _i; }

      LEDA_MEMORY(Face)
    };
  };



}; // HDS_items

template <typename R,class T> 
typename HDS_items::Vertex_wrapper<R,T>::Vertex::iv_iterator
HDS_items::Vertex_wrapper<R,T>::Vertex::nil_;

template <typename R,class T> 
typename HDS_items::Halfedge_wrapper<R,T>::Halfedge::fc_iterator
HDS_items::Halfedge_wrapper<R,T>::Halfedge::nil_;

#endif // CGAL_HDS_ITEMS_H

