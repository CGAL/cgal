// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_HDS_ITEMS_H
#define CGAL_HDS_ITEMS_H

#include <CGAL/license/Nef_2.h>


#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <list>
#include <boost/optional.hpp>
#include <boost/none.hpp>
#ifndef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <boost/any.hpp>
#endif

namespace CGAL {

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


//#ifndef CGAL_USE_LEDA
//#define LEDA_MEMORY(t) 
//#endif

template <typename Refs, typename Traits>
class Nef_vertex_2 {
public:
    typedef Refs           HalfedgeDS;
    typedef Nef_vertex_2   Base;
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
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    typedef void*  GenPtr;
    #else
    typedef boost::any GenPtr;
    #endif

    typedef typename Traits::Point Point;   // geometric embedding
    typedef typename Traits::Mark  Mark;    // information

    typedef typename std::list<Vertex_handle>::iterator iv_iterator;
private:
    Halfedge_handle              _h;
    Face_handle                  _f;
    Point                        _p;
    boost::optional<iv_iterator> _ivit;
    Mark                         _m;
    GenPtr                       _i;
public:

    Nef_vertex_2() : _h(),_f(),_ivit(),_m(0)
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      ,_i((GenPtr)0xABCD)
    #endif
    {}
    // constructs an uninitialized vertex concerning embedding 
    // and mark. All links are initialized by their default value.

    Nef_vertex_2( const Point& p) : 
        _h(),_f(),_p(p),_ivit(),_m(0)
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
        ,_i((GenPtr)0xABCD)
    #endif
    {}
    // constructs a vertex with embedding |p| and mark |m|.
    //  All links are initialized by their default value.

    bool is_isolated() const  { return _h == Halfedge_handle(); }
    // returns true iff |\Mvar| is isolated, else false.       
   
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

    iv_iterator ivit() const { return *_ivit; }
    void set_ivit(iv_iterator it) { _ivit = it; }
    void reset_ivit() { _ivit = boost::none; }
};


template <typename Refs, typename Traits>
class Nef_halfedge_2 {
public:
    typedef Refs           HalfedgeDS;
    typedef Nef_halfedge_2 Base;
    typedef Nef_halfedge_2 Base_base;
    typedef CGAL::Tag_true Supports_halfedge_prev;
    typedef CGAL::Tag_true Supports_halfedge_vertex;
    typedef CGAL::Tag_true Supports_halfedge_face;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Vertex                Vertex;
    typedef typename Refs::Face                  Face;
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    typedef void*  GenPtr;
    #else
    typedef boost::any GenPtr;
    #endif

    typedef typename std::list<Halfedge_handle>::iterator fc_iterator;

    typedef typename Traits::Mark  Mark;  // information 

protected:
    Halfedge_handle              opp, prv, nxt;
    Vertex_handle                _v;
    Face_handle                  _f;    
    boost::optional<fc_iterator> _fcit;
    Mark                         _m;
    GenPtr                       _i;
public:
       
    Nef_halfedge_2() : 
        opp(),prv(),nxt(),_v(),_f(),_fcit(),_m(0)
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
        ,_i((GenPtr)0xABCD)
    #endif
    {}
    /*{\Mcreate constructs an uninitialized halfedge concerning embedding 
      and mark. All links are initialized by their default value.}*/   

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

    fc_iterator fcit() const      { return *_fcit; }
    void set_fcit(fc_iterator it) { _fcit = it; }
    void reset_fcit()             { _fcit = boost::none; }

    bool is_hole_entry() const 
        /*{\Mop returns true iff |\Mvar| is entry point into a hole face
          cycle of |\Mvar.face()|.}*/
        { return !!_fcit; }
};

template <typename Refs, typename Traits>
class Nef_face_2 {
public:
    typedef Refs           HalfedgeDS;
    typedef Nef_face_2     Base;
    typedef CGAL::Tag_true Supports_face_halfedge;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
    typedef typename Refs::Vertex                Vertex;
    typedef typename Refs::Halfedge              Halfedge;
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    typedef void*  GenPtr;
    #else
    typedef boost::any GenPtr;
    #endif

    typedef typename Traits::Mark  Mark;  // mark information

    class Hole_iterator 
    /*{\Mtypemember iterator for face cycles. Fits the concept 
      |Halfedge_handle|.}*/ 
        : public std::list<Halfedge_handle>::iterator 
    {
        typedef typename std::list<Halfedge_handle>::iterator Ibase;
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
    {
        typedef typename std::list<Halfedge_handle>::const_iterator Ibase;
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
    {
        typedef typename std::list<Vertex_handle>::iterator Ibase;
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
    { 
        typedef typename std::list<Vertex_handle>::const_iterator Ibase;
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

    Nef_face_2() : _e(),_m(0)
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      ,_i((GenPtr)0xABCD) 
    #endif
    {}
    /*{\Mcreate constructs an uninitialized face with undefined mark,
      empty face cycle list, and empty isolated vertices list.}*/   

    ~Nef_face_2() { FC.clear(); IV.clear(); }

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
        FC.erase(h->fcit()); h->reset_fcit(); }

    void store_iv(Vertex_handle v)
        /*{\Mop stores vertex |v| as an isolated vertex of |\Mvar|.}*/
        { IV.push_back(v); v->set_ivit(--IV.end()); }

    void remove_iv(Vertex_handle v)
        /*{\Mop removes vertex |v| as an isolated vertex of |\Mvar|.
          \precond |v->is_isolated()| and |v| is stored in the
          isolated vertices list of |\Mvar|. 
          Postcondition: |!v->is_isolated()|.}*/
        { CGAL_assertion(v->is_isolated()); 
        IV.erase(v->ivit()); v->reset_ivit(); }
      
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
        {
        for (Hole_iterator hit = fc_begin(); hit!=fc_end(); ++hit) 
            hit->reset_fcit();
        for (Isolated_vertex_iterator vit = iv_begin(); vit!=iv_end(); ++vit) 
            vit->reset_ivit();
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
};

class HDS_items {
public:
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef Nef_vertex_2< Refs, Traits>    Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper {
        typedef Nef_halfedge_2< Refs, Traits>  Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef Nef_face_2< Refs, Traits>      Face;
    };
};

} // namespace CGAL

#endif // CGAL_HDS_ITEMS_H
