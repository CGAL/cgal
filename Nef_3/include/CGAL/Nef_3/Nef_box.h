// Copyright (c) 1997-2002,2005 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Peter Hachenberger    <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_NEF_BOX_H
#define CGAL_NEF_BOX_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/box_intersection_d.h>
#include <CGAL/Box_intersection_d/box_limits.h>

namespace CGAL {

template<class SNC_decorator>
class Nef_box : public Box_intersection_d::Box_d< double, 3 > {

  typedef typename SNC_decorator::SNC_structure          SNC_structure;
  typedef typename CGAL::SNC_intersection<SNC_structure> SNC_intersection;

  typedef typename SNC_decorator::Decorator_traits       Decorator_traits;
  typedef typename Decorator_traits::Vertex_handle        Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle      Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle     Halffacet_handle;
  typedef typename Decorator_traits::Halffacet_cycle_iterator
                                     Halffacet_cycle_iterator;
  typedef typename SNC_decorator::Infi_box                Infi_box;
  typedef typename SNC_decorator::Kernel                  Kernel;
  typedef typename SNC_decorator::Point_3                 Point_3;
  typedef typename Decorator_traits::SHalfedge_iterator   SHalfedge_iterator;
  typedef typename Decorator_traits::SHalfedge_around_facet_circulator
                                     SHalfedge_around_facet_circulator;

  typedef CGAL::SNC_const_decorator<SNC_structure>  Const_decorator;


  typedef std::pair<double, double> double_pair;
  typedef Box_intersection_d::box_limits<double> box_limits;

  Halffacet_handle f;
  Halfedge_handle  e;
  Vertex_handle v;
  enum Type { FACET, EDGE, VERTEX };
  Type type;

  void extend( const Point_3& p, const Tag_false& ) {
    std::pair<double, double> q[3];
    q[0] = CGAL::to_interval( p.x() );
    q[1] = CGAL::to_interval( p.y() );
    q[2] = CGAL::to_interval( p.z() );
    Box_intersection_d::Box_d< double, 3 >::extend(q);
  }

  void extend( const Point_3& p, const Tag_true& ) {
    double_pair q[3];
    if(Infi_box::degree(p.hx()) == 0)
      q[0] = CGAL::to_interval(p.x());
    else
      q[0] = p.x() > 0 ? double_pair(box_limits::sup(),box_limits::sup())
        : double_pair(box_limits::inf(),box_limits::inf());
    if(Infi_box::degree(p.hy()) == 0)
      q[1] = CGAL::to_interval(p.y());
    else
      q[1] = p.y() > 0 ? double_pair(box_limits::sup(),box_limits::sup())
        : double_pair(box_limits::inf(),box_limits::inf());
    if(Infi_box::degree(p.hz()) == 0)
      q[2] = CGAL::to_interval(p.z());
    else
      q[2] = p.z() > 0 ? double_pair(box_limits::sup(),box_limits::sup())
        : double_pair(box_limits::inf(),box_limits::inf());
    Box_intersection_d::Box_d< double, 3 >::extend(q);
  }

 public:
  Nef_box( Halffacet_handle f ) : f(f), type(FACET) {
    if( !Const_decorator::is_standard( f ) ) {
      init( true );
    } else {
      init( false );
#ifdef CGAL_NEF3_FACET_WITH_BOX
      std::pair<double, double> q[3];
      q[0] = CGAL::to_interval( f->b.min_coord(0) );
      q[1] = CGAL::to_interval( f->b.min_coord(1) );
      q[2] = CGAL::to_interval( f->b.min_coord(2) );
      Box_intersection_d::Box_d< double, 3 >::extend(q);
      q[0] = CGAL::to_interval( f->b.max_coord(0) );
      q[1] = CGAL::to_interval( f->b.max_coord(1) );
      q[2] = CGAL::to_interval( f->b.max_coord(2) );
      Box_intersection_d::Box_d< double, 3 >::extend(q);
#else
      Halffacet_cycle_iterator cycle_it = f->facet_cycles_begin();
      if( cycle_it.is_shalfedge() ) {
        SHalfedge_iterator edge_it(cycle_it);
        SHalfedge_around_facet_circulator
          start( edge_it ), end( edge_it );
        CGAL_For_all( start, end ) {
          const Point_3& p = start->prev()->source()->source()->point();
          extend( p, typename Is_extended_kernel<Kernel>::value_type());
        }
      } else
        CGAL_error_msg( "is facet first cycle a SHalfloop?");
#endif
    }
  }

  Nef_box( Halfedge_handle e ) :  e(e), type(EDGE) {

    if(!Const_decorator::is_standard(e->source() ) ||
       !Const_decorator::is_standard(e->twin()->source() ) ) {
      init( true );
    } else {
      init( false );
      extend( e->source()->point(),
              typename Is_extended_kernel<Kernel>::value_type());
      extend( e->twin()->source()->point(),
              typename Is_extended_kernel<Kernel>::value_type());
    }
  }

  Nef_box(Vertex_handle vin) : v(vin), type(VERTEX) {

    if(!Const_decorator::is_standard(v))
      init(true);
    else {
      init(false);
      extend(v->point(),
             typename Is_extended_kernel<Kernel>::value_type());
    }
  }

  Halffacet_handle get_halffacet() {
    CGAL_assertion( type == FACET );
    return f;
  }

  Halfedge_handle get_halfedge() {
    CGAL_assertion( type == EDGE );
    return e;
  }

  Vertex_handle get_vertex() {
    CGAL_assertion(type == VERTEX);
    return v;
  }
};

} //namespace CGAL

#endif // CGAL_NEF_BOX_H
