// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Meyer  <ameyer@mpi-sb.mpg.de>
#ifndef CGAL_BINOP_INTERSECTION_TESTS_H
#define CGAL_BINOP_INTERSECTION_TESTS_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_3/Nef_box.h>
#include <CGAL/Nef_3/Infimaximal_box.h>
#include <CGAL/Nef_3/SNC_const_decorator.h>
#include <vector>

namespace CGAL {

template<class SNC_decorator>
struct binop_intersection_test_segment_tree {
  typedef typename SNC_decorator::SNC_structure          SNC_structure;
  typedef typename CGAL::SNC_intersection<SNC_structure> SNC_intersection;

  typedef typename SNC_decorator::Decorator_traits        Decorator_traits;
  typedef typename Decorator_traits::Halfedge_iterator    Halfedge_iterator;
  typedef typename Decorator_traits::Halffacet_iterator   Halffacet_iterator;
  typedef typename SNC_decorator::Infi_box                Infi_box;
  typedef typename SNC_decorator::Point_3                 Point_3;

  typedef CGAL::SNC_const_decorator<SNC_structure>        Const_decorator;
  typedef CGAL::Nef_box<SNC_decorator>                    Nef_box;

  template<class Callback>
  struct Bop_edge0_face1_callback {
    Callback           &cb;

    Bop_edge0_face1_callback(Callback &cb)
    : cb(cb)
    {}

    void operator()( Nef_box& box0, Nef_box& box1 ) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif
      Halfedge_iterator  e0 = box0.get_halfedge();
      Halffacet_iterator f1 = box1.get_halffacet();
      if( Infi_box::degree( f1->plane().d() ) > 0 )
        return;
      Point_3 ip;
      if( SNC_intersection::does_intersect_internally( Const_decorator::segment(e0), f1, ip )) {
        cb(e0,f1,ip);
      }
    }
  };

  template<class Callback>
  struct Bop_edge1_face0_callback {
    Callback         &cb;

    Bop_edge1_face0_callback(Callback &cb)
    : cb(cb)
    {}

    void operator()( Nef_box& box0, Nef_box& box1 ) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif
      Halfedge_iterator  e1 = box0.get_halfedge();
      Halffacet_iterator f0 = box1.get_halffacet();
      if( Infi_box::degree( f0->plane().d() ) > 0 )
        return;
      Point_3 ip;
      if( SNC_intersection::does_intersect_internally( Const_decorator::segment( e1 ),
                                                       f0, ip ) )
        cb(e1,f0,ip);
    }
  };

  template<class Callback>
  struct Bop_edge0_edge1_callback  {
    Callback         &cb;

    Bop_edge0_edge1_callback(Callback &cb)
    : cb(cb)
    {}

    void operator()( Nef_box& box0, Nef_box& box1 ) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif
      Halfedge_iterator e0 = box0.get_halfedge();
      Halfedge_iterator e1 = box1.get_halfedge();
      Point_3 ip;
      if( SNC_intersection::does_intersect_internally( Const_decorator::segment( e0 ),
                                                       Const_decorator::segment( e1 ), ip ))
        cb(e0,e1,ip);
    }
  };

  template<class Callback>
  void operator()(Callback& cb0,
                  Callback& cb1,
                  const SNC_structure& snc0,
                  const SNC_structure& snc1)
  {
    Halfedge_iterator e0, e1;
    Halffacet_iterator f0, f1;
    std::vector<Nef_box> e0boxes, e1boxes, f0boxes, f1boxes;

    e0boxes.reserve(snc0.number_of_halfedges());
    e1boxes.reserve(snc1.number_of_halfedges());
    f0boxes.reserve(snc0.number_of_halffacets());
    f1boxes.reserve(snc1.number_of_halffacets());

    CGAL_forall_edges( e0, snc0) e0boxes.push_back( Nef_box( e0 ) );
    CGAL_forall_edges( e1, snc1) e1boxes.push_back( Nef_box( e1 ) );
    CGAL_forall_facets( f0, snc0) f0boxes.push_back( Nef_box( f0 ) );
    CGAL_forall_facets( f1, snc1) f1boxes.push_back( Nef_box( f1 ) );

    CGAL_NEF_TRACEN("start edge0 edge1");
    Bop_edge0_edge1_callback<Callback> callback_edge0_edge1( cb0 );
    box_intersection_d( e0boxes.begin(), e0boxes.end(),
                        e1boxes.begin(), e1boxes.end(),
                        callback_edge0_edge1);

    CGAL_NEF_TRACEN("start edge0 face1");
    Bop_edge0_face1_callback<Callback> callback_edge0_face1( cb0 );
    box_intersection_d( e0boxes.begin(), e0boxes.end(),
                        f1boxes.begin(), f1boxes.end(),
                        callback_edge0_face1);

    CGAL_NEF_TRACEN("start edge1 face0");
    Bop_edge1_face0_callback<Callback> callback_edge1_face0( cb1 );
    box_intersection_d( e1boxes.begin(), e1boxes.end(),
                        f0boxes.begin(), f0boxes.end(),
                        callback_edge1_face0);
  }
};

} //namespace CGAL

#endif
