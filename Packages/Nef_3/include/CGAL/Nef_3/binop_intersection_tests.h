// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Meyer  <ameyer@mpi-sb.mpg.de>
#ifndef CGAL_BINOP_INTERSECTION_TESTS_H
#define CGAL_BINOP_INTERSECTION_TESTS_H

#include <CGAL/box_intersection_d.h>
#include <vector>
#include <iostream>
#include <CGAL/Timer.h>

CGAL_BEGIN_NAMESPACE

template<class SNC_decorator>
struct binop_intersection_test_segment_tree {
  typedef typename SNC_decorator::SNC_structure          SNC_structure;
  typedef typename SNC_decorator::SNC_intersection       SNC_intersection;

  typedef typename SNC_decorator::Object_handle              Object_handle;
  typedef typename SNC_decorator::Vertex_const_iterator      Vertex_const_iterator;
  typedef typename SNC_decorator::Vertex_const_handle        Vertex_const_handle;
  typedef typename SNC_decorator::Halfedge_const_iterator    Halfedge_const_iterator;
  typedef typename SNC_decorator::Halfedge_const_handle      Halfedge_const_handle;
  typedef typename SNC_decorator::Halffacet_const_iterator   Halffacet_const_iterator;
  typedef typename SNC_decorator::Halffacet_const_handle     Halffacet_const_handle;
  typedef typename SNC_decorator::Halffacet_cycle_const_iterator
                                  Halffacet_cycle_const_iterator;
  typedef typename SNC_decorator::Infi_box               Infi_box;
  typedef typename SNC_decorator::Point_3                Point_3;
  typedef typename SNC_decorator::SHalfedge_const_handle     SHalfedge_const_handle;
  typedef typename SNC_decorator::SHalfedge_const_iterator   SHalfedge_const_iterator;
  typedef typename SNC_decorator::SHalfedge_around_facet_const_circulator
                                  SHalfedge_around_facet_const_circulator;

  typedef CGAL::SNC_const_decorator<SNC_structure>  Const_decorator;

  class Nef_box : public Box_intersection_d::Box_d< double, 3 >
  {
    Halffacet_const_handle f;
    Halfedge_const_handle  e;
    enum Type { FACET, EDGE };
    Type type;

    void extend( const Point_3& p ) {
      double q[3];
      q[0] = CGAL::to_double( p.x() );
      q[1] = CGAL::to_double( p.y() );
      q[2] = CGAL::to_double( p.z() );
      Box_intersection_d::Box_d< double, 3 >::extend(q);
    }

  public:
    Nef_box( Halffacet_const_handle f ) : f(f), type(FACET) {
      if( !Const_decorator::is_standard( f ) ) {
        init( true );
      } else {
        init( false );
        Halffacet_cycle_const_iterator cycle_it = f->facet_cycles_begin();
        if( cycle_it.is_shalfedge() ) {
	  SHalfedge_const_iterator edge_it(cycle_it);
          SHalfedge_around_facet_const_circulator
            start( edge_it ), end( edge_it );
          CGAL_For_all( start, end ) {
            const Point_3& p =
              SNC_decorator::point(
                SNC_decorator::source(
                  SNC_decorator::previous( start ) ) );
            extend( p );
          }
        } else
          CGAL_assertion_msg(0, "is facet first cycle a SHalfloop?");
      }
    }

    Nef_box( Halfedge_const_handle e ) :  e(e), type(EDGE)
    {
      if(!Const_decorator::is_standard( SNC_decorator::source( e ) ) ||
         !Const_decorator::is_standard( SNC_decorator::target( e ) ) )
      {
        init( true );
      } else {
        init( false );
        extend( SNC_decorator::point( SNC_decorator::source( e ) ) );
        extend( SNC_decorator::point( SNC_decorator::target( e ) ) );
      }
    }

    Halffacet_const_handle get_halffacet() {
      assert( type == FACET );
      return f;
    }

    Halfedge_const_handle get_halfedge() {
      assert( type == EDGE );
      return e;
    }
  };

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

    typedef std::pair<Halfedge_const_handle,
                      Halffacet_const_handle> Hash_key_type;


    Unique_hash_map< Hash_key_type,
                     bool,
                     Pair_hash_function > ignore;

    Bop_edge0_face1_callback(SNC_intersection &is, Callback &cb)
    : is(is), cb(cb), ignore( false )
    {}

    void operator()( Nef_box& box0, Nef_box& box1 ) {
      Halfedge_const_iterator  e0 = box0.get_halfedge();
      Halffacet_const_iterator f1 = box1.get_halffacet();
      if(ignore[ std::make_pair( e0, f1 ) ])
        return;
      if( Infi_box::degree( SNC_decorator::plane( f1 ).d() ) > 0 )
        return;
      Point_3 ip;
      if( is.does_intersect_internally( SNC_decorator::segment(e0), f1, ip )) {
        cb(e0,f1,ip);
        ignore[ std::make_pair( SNC_decorator::twin( e0 ), f1 ) ] = true;
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
      Halfedge_const_iterator  e1 = box0.get_halfedge();
      Halffacet_const_iterator f0 = box1.get_halffacet();
      if( Infi_box::degree( SNC_decorator::plane( f0 ).d() ) > 0 )
        return;
      Point_3 ip;
      if( is.does_intersect_internally( SNC_decorator::segment( e1 ),
                                        f0, ip ) )
        cb(e1,f0,ip);
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
      Halfedge_const_iterator e0 = box0.get_halfedge();
      Halfedge_const_iterator e1 = box1.get_halfedge();
      Point_3 ip;
      if( is.does_intersect_internally( SNC_decorator::segment( e0 ),
                                        SNC_decorator::segment( e1 ), ip ))
        cb(e0,e1,ip);
    }
  };

  template<class Callback>
  void operator()(Callback& cb0,
                  Callback& cb1,
                  const SNC_structure& sncp,
                  const SNC_structure& snc1i)
  {
    Halfedge_const_iterator e0, e1;
    Halffacet_const_iterator f0, f1;
    std::vector<Nef_box> a, b;
    SNC_intersection is( sncp );

    TRACEN("start edge0 edge1");
    Bop_edge0_edge1_callback<Callback> callback_edge0_edge1( is, cb0 );
    CGAL_forall_edges( e0, sncp)  a.push_back( Nef_box( e0 ) );
    CGAL_forall_edges( e1, snc1i) b.push_back( Nef_box( e1 ) );
    box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),
                        callback_edge0_edge1);
    a.clear();
    b.clear();

    TRACEN("start edge0 face1");
    Bop_edge0_face1_callback<Callback> callback_edge0_face1( is, cb0 );
    CGAL_forall_halfedges( e0, sncp ) a.push_back( Nef_box( e0 ) );
    CGAL_forall_facets( f1, snc1i)    b.push_back( Nef_box( f1 ) );
    box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),
                        callback_edge0_face1);
    a.clear();
    b.clear();

    TRACEN("start edge1 face0");
    Bop_edge1_face0_callback<Callback> callback_edge1_face0( is, cb1 );
    CGAL_forall_edges( e1, snc1i)  a.push_back( Nef_box( e1 ) );
    CGAL_forall_facets( f0, sncp ) b.push_back( Nef_box( f0 ) );
    box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),
                        callback_edge1_face0);
  }
};

CGAL_END_NAMESPACE

#endif
