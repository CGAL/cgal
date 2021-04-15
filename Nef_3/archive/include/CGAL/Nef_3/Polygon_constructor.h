// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF_POLYGON_CONSTRUCTOR_H
#define CGAL_NEF_POLYGON_CONSTRUCTOR_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Origin.h>
#include <iterator>
#include <iostream>

namespace CGAL {

template<class Nef3, typename forward_iterator>
class Polygon_constructor : public Modifier_base<typename Nef3::SNC_structure> {

  typedef Nef3                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_structure   SNC_structure;
  typedef typename SNC_structure::SM_decorator     SM_decorator;
  typedef typename SNC_structure::Vertex_handle    Vertex_handle;
  typedef typename SNC_structure::SVertex_handle   SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle     SFace_handle;
  typedef typename SNC_structure::Sphere_point     Sphere_point;
  typedef typename SNC_structure::Sphere_circle    Sphere_circle;
  typedef typename SNC_structure::Kernel           Kernel;
  typedef typename Kernel::Line_3                  Line_3;

  typedef typename std::iterator_traits<forward_iterator>::value_type
    point_iterator_pair;
  typedef typename point_iterator_pair::first_type
    point_iterator;

  forward_iterator begin, end;

 public:
  Polygon_constructor(forward_iterator begin_in, forward_iterator end_in) :
    begin(begin_in), end(end_in) {}

    void create_sphere_map(SNC_structure& snc,
                           point_iterator cur,
                           point_iterator prev,
                           point_iterator next,
                           Sphere_circle& c) {
      Vertex_handle v(snc.new_vertex(*cur, true));
      SM_decorator SM(&*v);
      SVertex_handle sv1(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
                                        true));
      SVertex_handle sv2(v->new_svertex(Sphere_point(ORIGIN+(*next-*cur)),
                                        true));
      SHalfedge_handle se(SM.new_shalfedge_pair(sv1,sv2));
      se->mark() = se->twin()->mark() = true;
      se->circle() = c;
      se->twin()->circle() = c.opposite();
      SFace_handle sf(v->new_sface());
      SM.link_as_face_cycle(se,sf);
    }

    Sphere_circle find_supporting_plane(point_iterator pbegin,
                                        point_iterator pend) {
      Line_3 l;
      point_iterator pnext(pbegin), pprev(pend), pcur(pbegin);
      --pprev;
      ++pnext;
      for(;pcur!=pend; ++pcur,++pprev,++pnext) {
        if(pprev == pend) pprev = pbegin;
        if(pnext == pend) pnext = pbegin;
        l = Line_3(*pcur,*pnext);
        if(!l.has_on(*pprev))
          return Sphere_circle(Sphere_point(ORIGIN+(*pprev-*pcur)),
                               Sphere_point(ORIGIN+(*pnext-*pcur)));
      }
      std::cerr << "all points lie on a common line" << std::endl;
      return Sphere_circle();
    }

 public:
    void operator()(SNC_structure& snc) {
      Sphere_circle c;
      point_iterator pbegin, pend, pnext, pprev;
      for(;begin != end; ++begin) {
        pprev = pend = begin->second;
        pnext = pbegin = begin->first;
        --pprev;
        ++pnext;
        c = find_supporting_plane(pbegin, pend);
        if(c == Sphere_circle()) continue;
        for(;pbegin!=pend; ++pbegin,++pprev, ++pnext) {
          if(pprev == pend) pprev = begin->first;
          if(pnext == pend) pnext = begin->first;
          create_sphere_map(snc,pbegin,pprev,pnext,c);
        }
      }
    }
};

} //namespace CGAL
#endif // CGAL_NEF_POLYGON_CONSTRUCTOR_H
