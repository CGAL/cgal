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
// Author(s)     : Andreas Meyer  <ameyer@mpi-sb.mpg.de>
#ifndef CGAL_BINOP_INTERSECTION_TESTS_H
#define CGAL_BINOP_INTERSECTION_TESTS_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_3/Nef_box.h>
#include <CGAL/Nef_3/Infimaximal_box.h>
#include <vector>
#include <iostream>
#include <CGAL/Timer.h>

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
    SNC_intersection   &is;
    Callback           &cb;

    struct Pair_hash_function {
      typedef std::size_t result_type;

      template <class H>
      std::size_t
      operator() (const H& h) const {
        return
          std::size_t(&*(h.first)) / sizeof
          (typename std::iterator_traits<typename H::first_type>::value_type)
              +
          std::size_t(&*(h.second)) / sizeof
          (typename std::iterator_traits<typename H::second_type>::value_type);
      }
    };

    Bop_edge0_face1_callback(SNC_intersection &is, Callback &cb)
    : is(is), cb(cb)
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
      if( is.does_intersect_internally( Const_decorator::segment(e0), f1, ip )) {
        cb(e0,make_object(f1),ip);
      }
    }
  };


  template<class Callback>
  struct Bop_edge1_face0_callback {
    SNC_intersection &is;
    Callback         &cb;

    Bop_edge1_face0_callback(SNC_intersection &is, Callback &cb)
    : is(is), cb(cb)
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
      if( is.does_intersect_internally( Const_decorator::segment( e1 ),
                                        f0, ip ) )
        cb(e1,make_object(f0),ip);
    }
  };

  template<class Callback>
  struct Bop_edge0_edge1_callback  {
    SNC_intersection &is;
    Callback         &cb;

    Bop_edge0_edge1_callback(SNC_intersection &is, Callback &cb)
    : is(is), cb(cb)
    {}

    void operator()( Nef_box& box0, Nef_box& box1 ) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif
      Halfedge_iterator e0 = box0.get_halfedge();
      Halfedge_iterator e1 = box1.get_halfedge();
      Point_3 ip;
      if( is.does_intersect_internally( Const_decorator::segment( e0 ),
                                        Const_decorator::segment( e1 ), ip ))
        cb(e0,make_object(e1),ip);
    }
  };

  template<class Callback>
  void operator()(Callback& cb0,
                  Callback& cb1,
                  SNC_structure& sncp,
                  SNC_structure& snc1i)
  {
    Halfedge_iterator e0, e1;
    Halffacet_iterator f0, f1;
    std::vector<Nef_box> a, b;
    SNC_intersection is( sncp );

    CGAL_NEF_TRACEN("start edge0 edge1");
    Bop_edge0_edge1_callback<Callback> callback_edge0_edge1( is, cb0 );
    CGAL_forall_edges( e0, sncp)  a.push_back( Nef_box( e0 ) );
    CGAL_forall_edges( e1, snc1i) b.push_back( Nef_box( e1 ) );
#ifdef CGAL_NEF3_BOX_INTERSECTION_CUTOFF
    box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),
                        callback_edge0_edge1,
			CGAL_NEF3_BOX_INTERSECTION_CUTOFF,);
#else
    box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),
                        callback_edge0_edge1);
#endif
    a.clear();
    b.clear();

    CGAL_NEF_TRACEN("start edge0 face1");
    Bop_edge0_face1_callback<Callback> callback_edge0_face1( is, cb0 );
    CGAL_forall_edges( e0, sncp ) a.push_back( Nef_box( e0 ) );
    CGAL_forall_facets( f1, snc1i)    b.push_back( Nef_box( f1 ) );
#ifdef CGAL_NEF3_BOX_INTERSECTION_CUTOFF
    box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),
                        callback_edge0_face1,
			CGAL_NEF3_BOX_INTERSECTION_CUTOFF);
#else
    box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),
                        callback_edge0_face1);
#endif
    a.clear();
    b.clear();

    CGAL_NEF_TRACEN("start edge1 face0");
    Bop_edge1_face0_callback<Callback> callback_edge1_face0( is, cb1 );
    CGAL_forall_edges( e1, snc1i)  a.push_back( Nef_box( e1 ) );
    CGAL_forall_facets( f0, sncp ) b.push_back( Nef_box( f0 ) );
#ifdef CGAL_NEF3_BOX_INTERSECTION_CUTOFF
    box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),
                        callback_edge1_face0,
			CGAL_NEF3_BOX_INTERSECTION_CUTOFF);
#else
    box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),
                        callback_edge1_face0);
#endif
  }
};

} //namespace CGAL

#endif
