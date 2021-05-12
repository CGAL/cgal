// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF_POLYLINE_CONSTRUCTOR_H
#define CGAL_NEF_POLYLINE_CONSTRUCTOR_H

#include <CGAL/Nef_3/SNC_indexed_items.h>

namespace CGAL {

template<typename Items, typename SNC_structure>
class Sphere_map_creator {
  typedef typename SNC_structure::SM_decorator     SM_decorator;
  typedef typename SNC_structure::Vertex_handle    Vertex_handle;
  typedef typename SNC_structure::SVertex_handle   SVertex_handle;
  typedef typename SNC_structure::SFace_handle     SFace_handle;
  typedef typename SNC_structure::Sphere_point     Sphere_point;

 public:
  Sphere_map_creator() {}

  template<typename point_iterator>
  void create_end_sphere_map(SNC_structure& snc,
                             point_iterator cur,
                             point_iterator prev) {
    Vertex_handle v(snc.new_vertex(*cur, true));
    SM_decorator SM(&*v);
    SVertex_handle sv(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
                                     true));
    SFace_handle sf(v->new_sface());
    SM.link_as_isolated_vertex(sv,sf);
  }

  template<typename point_iterator>
  void create_sphere_map(SNC_structure& snc,
                         point_iterator cur,
                         point_iterator prev,
                         point_iterator next) {
    Vertex_handle v(snc.new_vertex(*cur, true));
    SM_decorator SM(&*v);
    SVertex_handle sv1(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
                                      true));
    SVertex_handle sv2(v->new_svertex(Sphere_point(ORIGIN+(*next-*cur)),
                                      true));
    SFace_handle sf(v->new_sface());
    SM.link_as_isolated_vertex(sv1,sf);
    SM.link_as_isolated_vertex(sv2,sf);
    }
};

template<typename SNC_structure>
class Sphere_map_creator<CGAL::SNC_indexed_items, SNC_structure> {
  typedef typename SNC_structure::SM_decorator     SM_decorator;
  typedef typename SNC_structure::Vertex_handle    Vertex_handle;
  typedef typename SNC_structure::SVertex_handle   SVertex_handle;
  typedef typename SNC_structure::SFace_handle     SFace_handle;
  typedef typename SNC_structure::Sphere_point     Sphere_point;

  bool first;
  int index;
 public:
  Sphere_map_creator() : first(true) {}

    template<typename point_iterator>
    void create_end_sphere_map(SNC_structure& snc,
                               point_iterator cur,
                               point_iterator prev) {
      Vertex_handle v(snc.new_vertex(*cur, true));
      SM_decorator SM(&*v);
      SVertex_handle sv(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
                                       true));
      SFace_handle sf(v->new_sface());
      SM.link_as_isolated_vertex(sv,sf);
      if(first) {
        sv->set_index();
        index = sv->get_index();
        first = false;
      } else
        sv->set_index(index);
    }

    template<typename point_iterator>
    void create_sphere_map(SNC_structure& snc,
                           point_iterator cur,
                           point_iterator prev,
                           point_iterator next) {
      Vertex_handle v(snc.new_vertex(*cur, true));
      SM_decorator SM(&*v);
      SVertex_handle sv1(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
                                        true));
      SVertex_handle sv2(v->new_svertex(Sphere_point(ORIGIN+(*next-*cur)),
                                        true));
      SFace_handle sf(v->new_sface());
      SM.link_as_isolated_vertex(sv1,sf);
      SM.link_as_isolated_vertex(sv2,sf);
      sv1->set_index(index);
      sv2->set_index();
      index = sv2->get_index();
    }
};

template<class Nef3, typename forward_iterator>
class Polyline_constructor : public Modifier_base<typename Nef3::SNC_structure> {

  typedef Nef3                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_structure   SNC_structure;
  typedef typename SNC_structure::SM_decorator     SM_decorator;
  typedef typename SNC_structure::Vertex_handle    Vertex_handle;
  typedef typename SNC_structure::SVertex_handle   SVertex_handle;
  typedef typename SNC_structure::SFace_handle     SFace_handle;
  typedef typename SNC_structure::Sphere_point     Sphere_point;
  typedef typename SNC_structure::Items            Items;

  typedef typename std::iterator_traits<forward_iterator>::value_type
    point_iterator_pair;
  typedef typename point_iterator_pair::first_type
    point_iterator;

  forward_iterator begin, end;

 public:
  Polyline_constructor(forward_iterator begin_in, forward_iterator end_in) :
    begin(begin_in), end(end_in) {}

 public:
    void operator()(SNC_structure& snc) {
      point_iterator pbegin, pend, pnext, pprev;
      Sphere_map_creator<Items, SNC_structure> smc;
      for(;begin != end; ++begin) {
        pend = begin->second;
        pprev = pnext = pbegin = begin->first;
        ++pnext;
        CGAL_assertion(pnext != pend);
        smc.create_end_sphere_map(snc,pbegin,pnext);
        for(++pbegin,++pnext; pnext!=pend; ++pbegin,++pprev,++pnext)
          smc.create_sphere_map(snc,pbegin,pprev,pnext);
        smc.create_end_sphere_map(snc,pbegin,pprev);
      }
    }
};

} //namespace CGAL
#endif // CGAL_NEF_POLYLINE_CONSTRUCTOR_H
