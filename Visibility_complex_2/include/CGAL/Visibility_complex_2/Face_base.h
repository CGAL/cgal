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
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_FACE_BASE_H
#define CGAL_VISIBILITY_COMPLEX_2_FACE_BASE_H

#include <cmath>
#include <queue>
#include <iterator>

#include <CGAL/Visibility_complex_2/ccw_cw_traits.h>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {

// Kludge used to convert Border_const_iterator to Border_iterator
template<class Border_it1,class Border_it2> void copy_border_iterator
(Border_it1& a,const Border_it2& b) {
  a.f=b.f;
  a.e=b.e;
  a.beyond=b.beyond;
}


template< class Vc_ >
class Face_base {
public:
    typedef typename Vc_::Gt                   Gt;
    typedef typename Gt::Disk                  Disk;
    typedef typename Vc_::Disk_handle          Disk_handle;
    typedef typename Vc_::Vertex               Vertex;
    typedef typename Vc_::Vertex_handle        Vertex_handle;
    typedef typename Vc_::Edge                 Edge;
    typedef typename Vc_::Edge_handle          Edge_handle;
    typedef typename Vc_::Face                 Face;
    typedef typename Vc_::Face_handle          Face_handle;
    typedef typename Vc_::Vertex_const_handle  Vertex_const_handle;
    typedef typename Vc_::Edge_const_handle    Edge_const_handle;
    typedef typename Vc_::Face_const_handle    Face_const_handle;

private:
    Vertex_handle   sup_;
    Vertex_handle   inf_;

    Edge_handle back_view_;
    Edge_handle front_view_;

    Edge_handle     bottom_edge_;
    Edge_handle     top_edge_;

public:
    // CONSTRUCTORS  -----------------------------------------------------------
    Face_base() 
	: sup_(0)         , inf_(0),
	  back_view_(0)   , front_view_(0) ,
	  bottom_edge_(0) , top_edge_(0)
  { }

    Vertex_handle sup() { return sup_; }
    Vertex_handle inf() { return inf_; }
    Vertex_const_handle sup() const { return sup_; }
    Vertex_const_handle inf() const { return inf_; }
    void set_sup(const Vertex_handle& v);
    void set_inf(const Vertex_handle& v);
    // -------------------------------------------------------------------------
    Edge_const_handle front_view()  const    { return front_view_; }
    Edge_const_handle back_view()   const    { return back_view_;  }
    Edge_handle front_view()      { return front_view_; }
    Edge_handle back_view()       { return back_view_;  }
    void set_back_view (Edge_handle e) { back_view_  = e;    }
    void set_front_view(Edge_handle e) { front_view_ = e;    }

    void set_top_edge(const Edge_handle& e)    { top_edge_    = e;     }
    void set_bottom_edge(const Edge_handle& e) { bottom_edge_ = e;     }
    Edge_const_handle top_edge()      const          { return top_edge_;     }
    Edge_const_handle bottom_edge()   const          { return bottom_edge_;  }
    Edge_handle top_edge()                { return top_edge_;     }
    Edge_handle bottom_edge()             { return bottom_edge_;  }

private:
    template <class ep,class er,class fp,class from> 
    struct Border_iterator_wrapper {
      fp super;
      Border_iterator_wrapper(fp f) : super(f) {};
      class Border_iterator
        :public std::iterator<std::bidirectional_iterator_tag,
                              Edge,ptrdiff_t,ep,er> {
        template<class a,class b> friend 
          void copy_border_iterator(a&,const b&);
        fp f;
        ep e;
        bool beyond;
      public:
        Border_iterator(fp face,ep edge)
          :f(face),e(edge),beyond(true) {
          if (e) {
            beyond=(e->ul()!=face)&&((e->sign()?e->ur():e->dr())!=face)&&
              (e->dl()!=face);
          }
          if (beyond) e=0;
        };
        Border_iterator():f(0),e(0),beyond(true) {};

        Border_iterator(from x) {
          copy_border_iterator(*this,x);
        }

        bool operator==(const Border_iterator& a) {
          if (beyond) return a.beyond&&f==a.f;
          return (!a.beyond)&&f==a.f&&e==a.e;
        }
        
        bool operator!=(const Border_iterator& a) {
          return ! (operator==(a));
        }
        er       operator*()  const { 
          return *(operator->());
        }
        ep operator->() const {
          return e;
        }
        Border_iterator& operator++() { 
          if (!e) return *this;
          if (e->object()) {
            Vertex_handle v=e->sup();
            if (!v) 
              beyond=true;
            else {
              ep ee=0;
              do {
                if (!ee) ee=v->ccw_source_edge(); else ee=v->ccw_target_edge();
                fp fs[3]=
                  {ee->ul(),ee->sign()?ee->ur():ee->dr(),ee->dl()};
                for (int i=0;i<3;++i) {
                  if (f==fs[i]) {
                    e=ee;
                    return *this;
                  }
                }
              } while (ee!=v->ccw_target_edge());
              beyond=true;
            }
          } else beyond=true;
          return *this; }
        Border_iterator& operator--() { 
          if (!e) return *this;
          if (beyond)
            beyond=false;
          else {
            if (e->object()) {
              Vertex_handle v=e->inf();
              if (v) {
                ep ee=0;
                do {
                  if (!ee) ee=v->cw_source_edge(); else ee=v->cw_target_edge();
                  fp fs[3]=
                    {ee->ul(),ee->sign()?ee->ur():ee->dr(),ee->dl()};
                  for (int i=0;i<3;++i) {
                    if (f==fs[i]) {
                      e=ee;
                      return *this;
                    }
                  }
                } while (ee!=v->cw_target_edge());
              }
            }
          }
          return *this;}
        Border_iterator operator++(int) {
          Border_iterator s(*this);
          ++*this;
          return s;
        }
        Border_iterator operator--(int) {
          Border_iterator s(*this);
          --*this;
          return s;
        }
      }; 

      Border_iterator top_begin() const {
        Vertex_handle v=super->inf();
        if (v) {
          if (super==v->sup()) 
            return Border_iterator(super,
                                   v->ccw_source_edge());
          if (!v->is_constraint()) return Border_iterator(super,0);
          v=v->pi();
          if (super==v->target_cusp_face()) 
            return Border_iterator(super,
                                   v->target_cusp_edge());
          if (super==v->source_cusp_face()||v->source_cusp_face()==0) {
            if (v->is_left_xx())
              return Border_iterator(super,
                                     v->pi()->ccw_target_edge());
            else
              return Border_iterator(super,
                                     v->pi()->ccw_source_edge());            
          }
          CGAL_error();
        }
        return Border_iterator();
      }
      Border_iterator top_end() const {
        Vertex_handle v=super->sup();
        if (v) {
          if (super==v->inf())
            return ++Border_iterator(super,
                                     v->cw_target_edge());
          if (!v->is_constraint()) return Border_iterator(super,0);
          if (super==v->target_cusp_face()) 
            return ++Border_iterator(super,
                                     v->target_cusp_edge());
          if (super==v->source_cusp_face()||v->source_cusp_face()==0) {
            if (v->is_left_xx())
              return ++Border_iterator(super,
                                       v->cw_target_edge());
            else
              return ++Border_iterator(super,
                                       v->cw_source_edge()); 
          }
          CGAL_error();
        } 
        return Border_iterator(super,0);
      }
      Border_iterator bottom_begin() const {
        Vertex_handle v=super->inf();
        if (v) {
          if (super==v->sup())
            return Border_iterator(super,
                                   v->ccw_target_edge());
          if (!v->is_constraint()) return Border_iterator(super,0);
          v=v->pi();
          if (super==v->source_cusp_face()) 
            return Border_iterator(super,
                                   v->source_cusp_edge());
          if (super==v->target_cusp_face()||v->target_cusp_face()==0) {
            if (v->is_xx_left())
              return Border_iterator(super,
                                     v->pi()->ccw_target_edge());
            else
              return Border_iterator(super,
                                     v->pi()->ccw_source_edge());            
          }
          CGAL_error();
        } 
        return Border_iterator();
      }
      Border_iterator bottom_end() const {
        Vertex_handle v=super->sup();
        if (v) {
          if (super==v->inf()) 
            return ++Border_iterator(super,
                                     v->cw_source_edge());
          if (!v->is_constraint()) return Border_iterator(super,0);
          if (super==v->source_cusp_face()) 
            return ++Border_iterator(super,
                                     v->source_cusp_edge());
          if (super==v->target_cusp_face()||v->target_cusp_face()==0) {
            if (v->is_left_xx())
              return ++Border_iterator(super,
                                       v->cw_target_edge());
            else
              return ++Border_iterator(super,
                                       v->cw_source_edge());            
          }
          CGAL_error();
        } 
        return Border_iterator(super,0);
      }
    };
    struct Bogus_border_iterator {
      Face* f; Edge* e; bool beyond;
    };
    typedef Border_iterator_wrapper<Edge*,Edge&,Face*,Bogus_border_iterator> 
      Biw;
    typedef Border_iterator_wrapper<const Edge*,const Edge&,const Face*,
      typename Biw::Border_iterator> Biwc;

public:
    typedef typename Biwc::Border_iterator Border_const_iterator;
    typedef typename Biw::Border_iterator Border_iterator;
    Border_const_iterator top_begin() const {
      return Biwc(static_cast<const Face*>(this)).top_begin();
    }
    Border_const_iterator top_end() const {
      return Biwc(static_cast<const Face*>(this)).top_end();
    }
    Border_const_iterator bottom_begin() const {
      return Biwc(static_cast<const Face*>(this)).bottom_begin();
    }
    Border_const_iterator bottom_end() const {
      return Biwc(static_cast<const Face*>(this)).bottom_end();
    }
    Border_iterator top_begin() {
      return Biw(static_cast<Face*>(this)).top_begin();
    }
    Border_iterator top_end() {
      return Biw(static_cast<Face*>(this)).top_end();
    }
    Border_iterator bottom_begin() {
      return Biw(static_cast<Face*>(this)).bottom_begin();
    }
    Border_iterator bottom_end() {
      return Biw(static_cast<Face*>(this)).bottom_end();
    }

    ~Face_base() {
      Edge_handle constraint_edge=0;
      Border_iterator border_end(static_cast<Face*>(this),0);
      if (inf()) {
        for (Border_iterator i=bottom_begin();
             i!=border_end;
             ++i) {
          // if i is a constraint edge, it must not be cleared of its face
          // yet, in order to keep the (top|bottom)_(begin|end) working.
          if (!i->object()) {
            constraint_edge=i.operator->();
            break;
          }
          if (i->dl()==static_cast<Face*>(this))
            i->set_adjacent_faces(0,i->sign()?i->ur():i->dr(),i->ul());
          else if ((i->sign()?i->ur():i->dr())==static_cast<Face*>(this))
            i->set_adjacent_faces(i->dl(),0,i->ul());
          else if (i->ul()==static_cast<Face*>(this))
            i->set_adjacent_faces(i->dl(),i->sign()?i->ur():i->dr(),0);
        }
        for (Border_iterator i=top_begin();
             i!=border_end;
             ++i) {
          if (!i->object()) {
            constraint_edge=i.operator->();
            break;
          }
          if (i->dl()==static_cast<Face*>(this))
            i->set_adjacent_faces(0,i->sign()?i->ur():i->dr(),i->ul());
          else if ((i->sign()?i->ur():i->dr())==static_cast<Face*>(this))
            i->set_adjacent_faces(i->dl(),0,i->ul());
          else if (i->ul()==static_cast<Face*>(this))
            i->set_adjacent_faces(i->dl(),i->sign()?i->ur():i->dr(),0);
        }     
      }
      if (sup()) {
//         if (!top_edge()) set_top_edge(sup()->cw_target_edge());
//         if (!bottom_edge()) set_bottom_edge(sup()->cw_source_edge());
        Border_iterator top=top_end();
        Border_iterator bot=bottom_end();
        if (top.operator->()) set_top_edge(top.operator->());
        if (bot.operator->()) set_bottom_edge(bot.operator->());
      }
      if (top_edge()){
        Edge_handle prev=0;
        for (Border_iterator i=
               Border_iterator(static_cast<Face*>(this),top_edge());
             i!=top_end()&&i.operator->()!=prev;--i) {
          if (!i->object()) {
            constraint_edge=i.operator->();
            break;
          }
          if (i->dl()==static_cast<Face*>(this))
            i->set_adjacent_faces(0,i->sign()?i->ur():i->dr(),i->ul());
          else if ((i->sign()?i->ur():i->dr())==static_cast<Face*>(this))
            i->set_adjacent_faces(i->dl(),0,i->ul());
          else if (i->ul()==static_cast<Face*>(this))
            i->set_adjacent_faces(i->dl(),i->sign()?i->ur():i->dr(),0);
          prev=i.operator->();
        }
      }
      if (bottom_edge()) {
        Edge_handle prev=0;
        for (Border_iterator i=
               Border_iterator(static_cast<Face*>(this),bottom_edge());
             i!=bottom_end()&&i.operator->()!=prev;--i) {
          if (!i->object()) {
            constraint_edge=i.operator->();
            break;
          }
          if (i->dl()==static_cast<Face*>(this))
            i->set_adjacent_faces(0,i->sign()?i->ur():i->dr(),i->ul());
          else if ((i->sign()?i->ur():i->dr())==static_cast<Face*>(this))
            i->set_adjacent_faces(i->dl(),0,i->ul());
          else if (i->ul()==static_cast<Face*>(this))
            i->set_adjacent_faces(i->dl(),i->sign()?i->ur():i->dr(),0);
          prev=i.operator->();
        }
      }
      if (constraint_edge) {
        constraint_edge->set_adjacent_faces(0,0,0);
      }
      if (sup() != 0 && static_cast<Face*>(this) == sup()->inf()) 
        sup()->set_inf(0);	
      if (inf() != 0 && static_cast<Face*>(this) == inf()->sup())
        inf()->set_sup(0);	
    }
};

// -----------------------------------------------------------------------------

template < class Vc_ >
void 
Face_base<Vc_>::set_sup(const Vertex_handle& v)
{
    sup_ = v;
    if (v != 0 && !v->is_constraint()) 
	v->set_inf(static_cast<Face_handle>(this));
}

template < class Vc_ >
void 
Face_base<Vc_>::set_inf(const Vertex_handle& v)
{
    inf_ = v;
    if (v != 0 && !v->is_constraint()) 
	v->set_sup(static_cast<Face_handle>(this));
}

// -----------------------------------------------------------------------------
}
CGAL_END_NAMESPACE

#endif // CGAL_VISIBILITY_COMPLEX_2_FACE_BASE_H
