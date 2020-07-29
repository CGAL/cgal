// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de>
//                 Miguel Granados    <granados@mpi-sb.mpg.de>
//                 Susan Hert         <hert@mpi-sb.mpg.de>
//                 Lutz Kettner       <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenb@mpi-sb.mpg.de>
#ifndef CGAL_NEF3_BINARY_OPERATION_H
#define CGAL_NEF3_BINARY_OPERATION_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Nef_S2/Normalizing.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_3/SNC_external_structure.h>
#include <CGAL/Nef_3/SNC_point_locator.h>
#include <CGAL/Nef_3/binop_intersection_tests.h>
#include <CGAL/Nef_3/ID_support_handler.h>
//#include <CGAL/Nef_3/Edge_edge_overlay.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 19
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

#if defined (CGAL_NEF3_TIMER_OVERLAY) || defined(CGAL_NEF3_TIMER_INTERSECTION)
CGAL::Timer timer_overlay;
#endif

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
int number_of_point_location_queries;
int number_of_ray_shooting_queries;
CGAL::Timer timer_point_location;
CGAL::Timer timer_ray_shooting;
#endif

#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
int number_of_edge_facet_overlays=0;
int number_of_clones=0;
int number_of_sphere_sweeps=0;
CGAL::Timer timer_sphere_sweeps;
#endif

#ifdef CGAL_NEF3_TIMER_PLANE_SWEEPS
int number_of_plane_sweeps=0;
CGAL::Timer timer_plane_sweeps;
#endif

#ifdef CGAL_NEF3_DUMP_STATISTICS
int number_of_intersections;
int number_of_intersection_candidates;
#endif

template <typename Map>
class Binary_operation : public CGAL::SNC_decorator<Map> {
 public:
  typedef Map SNC_structure;
  typedef typename SNC_structure::Items                Items;
  typedef typename Map::Sphere_map                     Sphere_map;
  typedef CGAL::SNC_decorator<SNC_structure>           SNC_decorator;
  typedef SNC_decorator                                Base;
  typedef CGAL::SNC_constructor<Items, SNC_structure>  SNC_constructor;
  typedef CGAL::SNC_external_structure<Items, SNC_structure>
    SNC_external_structure;
  typedef CGAL::SM_decorator<Sphere_map>               SM_decorator;
  typedef CGAL::SNC_SM_overlayer<Items, SM_decorator>  SM_overlayer;
  typedef CGAL::SM_point_locator<SM_decorator>         SM_point_locator;
  typedef CGAL::SNC_point_locator<SNC_decorator>       SNC_point_locator;

  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Volume_handle Volume_handle;
  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;
  typedef typename SNC_structure::SFace_handle SFace_handle;

  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;
  typedef typename SNC_structure::Volume_iterator Volume_iterator;
  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;
  typedef typename SNC_structure::SFace_iterator SFace_iterator;

  typedef typename SNC_structure::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Shell_entry_iterator Shell_entry_iterator;

  typedef typename SNC_structure::Object_handle Object_handle;

  typedef typename Base::Vertex_const_handle Vertex_const_handle;
  typedef typename Base::Volume_const_handle Volume_const_handle;

  typedef typename Base::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Base::SHalfedge_const_iterator SHalfedge_const_iterator;
  typedef typename Base::SHalfloop_const_iterator SHalfloop_const_iterator;

  typedef typename Base::Point_3 Point_3;
  typedef typename Base::Plane_3 Plane_3;

  typedef typename Base::Mark Mark;

  typedef CGAL::ID_support_handler<Items, SNC_decorator> Association;

 public:
  Binary_operation(SNC_structure& W) : Base(W) {}

  template <typename Selection, typename Association>
  Vertex_handle binop_local_views( Vertex_const_handle v0, Vertex_const_handle v1,
                                   const Selection& BOP, SNC_structure& rsnc
                                   ,Association& A)
    /*{\opOverlays two spheres maps.}*/ {
    //    CGAL_NEF_SETDTHREAD(19*43*131);

    CGAL_assertion( v0->point() == v1->point());
    Vertex_handle v01 = rsnc.new_vertex( v0->point(), BOP( v0->mark(),v1->mark()));
    //    std::cerr <<"BOP Vertex "<< v0->point() << ":"
    //                     << v0->mark()<<" "<<v1->mark()<<std::endl;
    CGAL_NEF_TRACEN("  binop result on vertex "<<&*v01<<" on "<<&*(v01->sncp()));
    SM_overlayer O(&*v01);
    O.subdivide( &*v0, &*v1, A);
    O.select( BOP);
    O.simplify(A);

    return v01;
  }

  Vertex_handle create_local_view_on( const Point_3& p, Halfedge_handle e) {
    SNC_constructor C(*this->sncp());
    return C.create_from_edge( e, p);
  }

  Vertex_handle create_local_view_on( const Point_3& p, Halffacet_handle f) {
    SNC_constructor C(*this->sncp());
    return C.create_from_facet( f, p);
  }

  Vertex_handle create_local_view_on( const Point_3& p, Volume_const_handle c) {
    Vertex_handle v = this->sncp()->new_vertex( p, c->mark());
    SM_decorator SD(&*v);
    SFace_handle f = SD.new_sface();
    f->mark() = c->mark();
    CGAL_NEF_TRACEN("volume "<<&*c<<" marked as "<<c->mark());
    return v;
  }

  template <typename SNC_decorator,
            typename Selection,
            typename Association>
  class Intersection_call_back :
    public SNC_point_locator::Intersection_call_back
  {
    typedef typename SNC_decorator::Decorator_traits Decorator_traits;
    typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
    typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;

  public:
    Intersection_call_back( SNC_structure& s0, SNC_structure& s1,
                            const Selection& _bop, SNC_structure& r,
                            bool invert_order, Association& Ain) :
      snc0(s0), snc1(s1), bop(_bop), result(r),
      inverse_order(invert_order), A(Ain) {}

      void operator()(Halfedge_handle e0, Object_handle o1, const Point_3& ip)
      const {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersections;
#endif

      Halfedge_handle e;
      Halffacet_handle f;

      Point_3 p(normalized(ip));

      CGAL_NEF_TRACEN("Intersection_call_back: intersection reported on " << p << " (normalized: " << normalized(p) << " )");
#ifdef CGAL_NEF_DEBUG
      CGAL_NEF_TRACEN("edge 0 has source " << e0->source()->point() << " and direction " << e0->vector());
      if( CGAL::assign( e, o1)) {
        CGAL_NEF_TRACEN("edge 1 has source " << e->source()->point() << " and direction " << e->vector());
      }
      else if( CGAL::assign( f, o1)) {
        CGAL_NEF_TRACEN("face 1 has plane equation " << f->plane());
      }
      else
              CGAL_error_msg( "wrong handle");
#endif

#if defined (CGAL_NEF3_TIMER_OVERLAY) || (CGAL_NEF3_TIMER_INTERSECTION)
      timer_overlay.start();
#endif

      if( CGAL::assign( e, o1)) {
        //        std::cerr << "inverse order " << inverse_order << std::endl;

#ifdef CGAL_NEF_EXPERIMENTAL_CODE
        typename CGAL::Edge_edge_overlay<SNC_structure> eeo(result, e0, e);
        Sphere_map* M0 = eeo.create_edge_edge_overlay(p, bop, inverse_order, A);
        SM_overlayer O(M0);
        O.simplify(A);
#else
        Binary_operation D(result);
        Vertex_handle v0, v1;
        v0 = D.create_local_view_on( p, e0);
        v1 = D.create_local_view_on( p, e);
        if( inverse_order)
          std::swap( v0, v1);
        D.binop_local_views( v0, v1, bop, result,A);
        result.delete_vertex(v0);
        result.delete_vertex(v1);
#endif
      }
      else if( CGAL::assign( f, o1)) {
#ifdef CGAL_NEF3_OVERLAY_BY_HAND_OFF
        Binary_operation D(result);
        Vertex_handle v0, v1;
        v0 = D.create_local_view_on( p, e0);
        v1 = D.create_local_view_on( p, f);
        if( inverse_order)
          std::swap( v0, v1);
        D.binop_local_views( v0, v1, bop, result,A);
        result.delete_vertex(v0);
        result.delete_vertex(v1);
#else // CGAL_NEF3_OVERLAY_BY_HAND_OFF
        SNC_constructor C(result);
        Sphere_map* M0 = C.create_edge_facet_overlay(e0, f, p, bop, inverse_order, A);
        SM_overlayer O(M0);
        O.simplify(A);
#endif // CGAL_NEF3_OVERLAY_BY_HAND_OFF
      }
      else
        CGAL_error_msg( "wrong handle");

#if defined (CGAL_NEF3_TIMER_OVERLAY) || (CGAL_NEF3_TIMER_INTERSECTION)
      timer_overlay.stop();
#endif

    }
  private:
    const SNC_structure& snc0;
    const SNC_structure& snc1;
    const Selection& bop;
    SNC_structure& result;
    bool inverse_order;
    Association& A;
  };

  template <typename Selection>
    void operator()( SNC_point_locator* pl0,
                     const SNC_structure& snc1,
                     const SNC_point_locator* pl1,
                     const SNC_structure& snc2,
                     const SNC_point_locator* pl2,
                     const Selection& BOP)
      /*{\opPerforms a binary operation defined on |BOP| between two
      SNC structures.  The input structures are not modified and the
      result of the operation is stored in |result|.
      \precondition: the structure |result| is empty.}*/
  {
    //    CGAL_NEF_SETDTHREAD(23);
    CGAL_assertion( this->sncp()->is_empty());
    CGAL_assertion( pl1 != nullptr && pl2 != nullptr);
    //    CGAL_NEF_SETDTHREAD(19*13*43*37);

#ifdef CGAL_NEF3_TIMER_BINARY_OPERATION
    CGAL::Timer timer_binary_operation;
    timer_binary_operation.start();
#endif

#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
    timer_sphere_sweeps.reset();
    number_of_sphere_sweeps=0;
    number_of_edge_facet_overlays=0;
    number_of_clones=0;
#endif

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
    timer_point_location.reset();
    number_of_point_location_queries=0;
#endif

#ifdef CGAL_NEF3_TIMER_OVERLAY
      timer_overlay.reset();
#endif

#ifdef CGAL_NEF3_DUMP_STATISTICS
      number_of_intersections=0;
      number_of_intersection_candidates=0;
#endif

    Unique_hash_map<Vertex_const_handle, bool> ignore(false);
    Vertex_const_iterator v0;

    //    CGAL_NEF_SETDTHREAD(19*43*131);
    CGAL_NEF_TRACEN("=> binary operation");

#ifdef CGAL_NEF3_FACET_WITH_BOX
    SNC_constructor C1(snc1);
    C1.create_box();
    SNC_constructor C2(snc2);
    C2.create_box();
#endif

    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "
                    << this->sncp()->number_of_vertices());

    CGAL_NEF_TRACEN("=> for all v0 in snc1, qualify v0 with respect snc2");
    //    int i=2;
    Association A;
    SHalfedge_const_iterator sei;
    CGAL_forall_shalfedges(sei, snc1)
      A.initialize_hash(sei);
    CGAL_forall_shalfedges(sei, snc2)
      A.initialize_hash(sei);
    SHalfloop_const_iterator sli;
    CGAL_forall_shalfloops(sli, snc1)
      A.initialize_hash(sli);
    CGAL_forall_shalfloops(sli, snc2)
      A.initialize_hash(sli);

    CGAL_forall_vertices( v0, snc1) {
      CGAL_assertion(!ignore[v0]);
      Point_3 p0(v0->point());
      Vertex_handle v;
      Halfedge_handle e;
      Halffacet_handle f;
      Volume_handle c;
      CGAL_NEF_TRACEN("Locating point " << p0);

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
      ++number_of_point_location_queries;
      timer_point_location.start();
#endif
      Object_handle o = pl2->locate(p0);
#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
      timer_point_location.stop();
#endif

#if defined(CGAL_NEF3_TIMER_OVERLAY)
      timer_overlay.start();
#endif
      if( CGAL::assign( v, o)) {
        CGAL_NEF_TRACEN("p0 found on vertex");
        binop_local_views( v0, v, BOP, *this->sncp(),A);
        ignore[v] = true;
      }
      else if( CGAL::assign( e, o)) {
        CGAL_NEF_TRACEN("p0 found on edge");
        Vertex_handle v1 = create_local_view_on( p0, e);
        binop_local_views( v0, v1, BOP, *this->sncp(),A);
        this->sncp()->delete_vertex(v1);
      }
      else if( CGAL::assign( f, o)) {
        CGAL_NEF_TRACEN("p0 found on facet" << f->plane());
        Vertex_handle v1 = create_local_view_on( p0, f);
        binop_local_views( v0, v1, BOP, *this->sncp(),A);
        this->sncp()->delete_vertex(v1);
      }
      else if( CGAL::assign( c, o)) {
        CGAL_NEF_TRACEN("p0 found on volume with mark " << c->mark());
#ifdef CGAL_NEF3_OVERLAY_IF_NEEDED_OFF
        if(true) {
#else
        if( BOP( true, c->mark()) != BOP( false, c->mark())) {
#endif
#ifdef CGAL_NEF3_OVERLAY_BY_HAND_OFF
          Vertex_handle v1 = create_local_view_on( p0, c);
          binop_local_views( v0, v1, BOP, *this->sncp(),A);
          this->sncp()->delete_vertex(v1);
#else
          SNC_constructor C(*this->sncp());
          Vertex_handle v1 = C.clone_SM(v0);
          SM_decorator SM(&*v1);
          SM.change_marks(BOP, c->mark());
          SM_overlayer O(&*v1);
          O.simplify(A);
#endif
        } else {
          CGAL_NEF_TRACEN("vertex in volume deleted " << std::endl <<
                 "  vertex: " <<  v0->point() << std::endl <<
                 "  mark of volume: " << c->mark());
        }
      }
      else CGAL_error_msg( "wrong handle");

#if defined(CGAL_NEF3_TIMER_OVERLAY)
      timer_overlay.stop();
#endif
      }
    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "
                    << this->sncp()->number_of_vertices());

    CGAL_NEF_TRACEN("=> for all v1 in snc1, qualify v1 with respect snc0");
    CGAL_forall_vertices( v0, snc2) {

      if(ignore[v0]) continue;
      Point_3 p1(v0->point());
      Halfedge_handle e;
      Halffacet_handle f;
      Volume_handle c;
      CGAL_NEF_TRACEN("Locating point " << p1);

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
      number_of_point_location_queries++;
      timer_point_location.start();
#endif
      Object_handle o = pl1->locate(p1);
#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
      timer_point_location.stop();
#endif

      CGAL_assertion_code(Vertex_handle v);
      CGAL_assertion( !CGAL::assign( v, o));

#if defined(CGAL_NEF3_TIMER_OVERLAY)
      timer_overlay.start();
#endif
      if( CGAL::assign( e, o)) {
        CGAL_NEF_TRACEN("p1 found on edge");
        Vertex_handle v1 = create_local_view_on( p1, e);
        binop_local_views( v1, v0, BOP, *this->sncp(),A);
        this->sncp()->delete_vertex(v1);
      }
      else if( CGAL::assign( f, o)) {
        CGAL_NEF_TRACEN("p1 found on facet");
        Vertex_handle v1 = create_local_view_on( p1, f);
        binop_local_views( v1, v0, BOP, *this->sncp(),A);
        this->sncp()->delete_vertex(v1);
      }
      else if( CGAL::assign( c, o)) {
        CGAL_NEF_TRACEN("p1 found on volume with mark " << c->mark());
#ifdef CGAL_NEF3_OVERLAY_IF_NEEDED_OFF
        if(true)
#else
        if( BOP( c->mark(), true) != BOP( c->mark(), false))
#endif
        {
#ifdef CGAL_NEF3_OVERLAY_BY_HAND_OFF
          Vertex_handle v1 = create_local_view_on( p1, c);
          binop_local_views( v1, v0, BOP, *this->sncp(),A);
          this->sncp()->delete_vertex(v1);
#else
          SNC_constructor C(*this->sncp());
          Vertex_handle v1 = C.clone_SM(v0);
          SM_decorator SM(&*v1);
          SM.change_marks(c->mark(), BOP);
          SM_overlayer O(&*v1);
          O.simplify(A);
#endif
        } else {
          CGAL_NEF_TRACEN("vertex in volume deleted " << std::endl <<
                          "  vertex: " <<  v0->point() << std::endl <<
                          "  mark of volume: " << c->mark());
        }
      }
      else CGAL_error_msg( "wrong handle");

#if defined(CGAL_NEF3_TIMER_OVERLAY)
      timer_overlay.stop();
#endif
    }

    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "<<
                    this->sncp()->number_of_vertices());

    // Each time the intersect method of the point locator for the
    // SNC structure finds an intersection between the segment defined
    // by an edge on the other SNC structure, the call back method is
    // called with the intersecting objects and the intersection point.
    // The responsability of the call back functor is to construct the
    // local view on the intersection point on both SNC structures,
    // overlay them and add the resulting sphere map to the result.

    // CGAL_NEF_SETDTHREAD(19*509*43*131);

    Intersection_call_back<SNC_decorator, Selection, Association> call_back0
      ( const_cast<SNC_structure&>(snc1), const_cast<SNC_structure&>(snc2),
        BOP, *this->sncp(), false, A);
    Intersection_call_back<SNC_decorator, Selection, Association> call_back1
      ( const_cast<SNC_structure&>(snc2), const_cast<SNC_structure&>(snc2),
        BOP, *this->sncp(), true, A);

#ifdef CGAL_NEF3_TIMER_INTERSECTION
    double split_intersection = timer_overlay.time();
    CGAL::Timer timer_intersection;
    timer_intersection.start();
#endif

    // choose between intersection algorithms
#ifdef CGAL_NEF3_INTERSECTION_BY_KDTREE
    Halfedge_iterator e0, e1;
    /*
    CGAL_NEF_TRACEN("=> finding edge-edge intersections...");
    CGAL_forall_edges( e0, snc1) {
      //      ee_intersections++;
      pl2->intersect_with_edges( e0, call_back0);
    }

    CGAL_NEF_TRACEN("number of vertices (so far...) = "
    << this->sncp()->number_of_vertices());
    CGAL_NEF_TRACEN("=> finding edge0-facet1 intersections...");
    CGAL_forall_edges( e0, snc1) {
      //      ef_intersections++;
      pl2->intersect_with_facets( e0, call_back0);
    }
    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "
    << this->sncp()->number_of_vertices());
    */

    CGAL_forall_edges(e0,const_cast<SNC_structure&>(snc1))
      pl2->intersect_with_edges_and_facets(e0,call_back0);


    CGAL_NEF_TRACEN("=> finding edge1-facet0 intersections...");
    CGAL_forall_edges( e1,const_cast<SNC_structure&>(snc2)) {
      //      ef_intersections++;
      pl1->intersect_with_facets( e1, call_back1);
    }
    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "
                    << this->sncp()->number_of_vertices());
#elif defined CGAL_NEF3_INTERSECTION_NAIVE

    CGAL::SNC_point_locator_naive<SNC_decorator> pln1;
    CGAL::SNC_point_locator_naive<SNC_decorator> pln2;
    pln1.initialize(const_cast<SNC_structure*>(&snc1));
    pln2.initialize(const_cast<SNC_structure*>(&snc2));
    Halfedge_iterator e0;
    CGAL_forall_edges(e0,const_cast<SNC_structure&>(snc1))
      pln2.intersect_with_edges_and_facets(e0,call_back0);
    CGAL_forall_edges(e0,const_cast<SNC_structure&>(snc2))
      pln1.intersect_with_facets( e0, call_back1);

#else
    CGAL_NEF_TRACEN("intersection by fast box intersection");
        binop_intersection_test_segment_tree<SNC_decorator> binop_box_intersection;
        binop_box_intersection(call_back0, call_back1,
                               const_cast<SNC_structure&>(snc1),
                               const_cast<SNC_structure&>(snc2));
#endif

#ifdef CGAL_NEF3_TIMER_INTERSECTION
    timer_intersection.stop();
    if(cgal_nef3_timer_on)
      std::cout << "Runtime_intersection: "
                << timer_intersection.time()-timer_overlay.time()+split_intersection
                << std::endl;
#endif

    CGAL_NEF_TRACEN("=> resultant vertices (before simplification): ");
    CGAL_assertion_code(CGAL_forall_vertices( v0, *this->sncp())
                          CGAL_NEF_TRACEN(&*v0<<" "<<v0->point()));

    SNC_external_structure es(*this->sncp(), pl0);
    es.build_after_binary_operation(A);

#ifdef CGAL_NEF3_TIMER_PLANE_SWEEPS
    if(cgal_nef3_timer_on) {
      std::cout << "Number_of_plane_sweeps: "
                << number_of_plane_sweeps << std::endl;
      std::cout << "Runtime_plane_sweeps: "
                << timer_plane_sweeps.time() << std::endl;
    }
#endif

#ifdef CGAL_NEF3_DUMP_STATISTICS
    if(cgal_nef3_timer_on) {
      std::cout << "Vertices_in_object_A: "
                << snc1.number_of_vertices() << std::endl;
      std::cout << "Vertices_in_object_B: "
                << snc2.number_of_vertices() << std::endl;
      std::cout << "Number_of_intersections: "
                << number_of_intersections << std::endl;
      std::cout << "Number_of_intersection_candidates: "
                << number_of_intersection_candidates << std::endl;
      std::cout << "Vertices_in_Result: "
                << this->sncp()->number_of_vertices() << std::endl;
    }
#endif

#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
    if(cgal_nef3_timer_on) {
      std::cout << "Number_of_edge_facet_overlays: "
                << number_of_edge_facet_overlays << std::endl;
      std::cout << "Number_of_clones: "
                << number_of_clones << std::endl;
      std::cout << "Number_of_sphere_sweeps: "
                << number_of_sphere_sweeps << std::endl;
      std::cout << "Runtime_sphere_sweeps: "
                << timer_sphere_sweeps.time() << std::endl;
    }
#endif

#if defined (CGAL_NEF3_TIMER_SPHERE_SWEEPS) && defined (CGAL_NEF3_TIMER_PLANE_SWEEPS)
    if(cgal_nef3_timer_on) {
      std::cout << "Runtime_all_sweeps: "
                << timer_sphere_sweeps.time()+timer_plane_sweeps.time() << std::endl;
    }
#endif

#ifdef CGAL_NEF3_TIMER_OVERLAY
    if(cgal_nef3_timer_on) {
      std::cout << "Runtime_overlay: "
                << timer_overlay.time() << std::endl;
    }
#endif

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
    if(cgal_nef3_timer_on) {
      std::cout << "Number_of_ray_shooting_queries: "
                << number_of_ray_shooting_queries << std::endl;
      std::cout << "Runtime_ray_shooting: "
                << timer_ray_shooting.time() << std::endl;
      std::cout << "Number_of_point_location_queries: "
                << number_of_point_location_queries << std::endl;
      std::cout << "Runtime_point_location: "
                << timer_point_location.time() << std::endl;
      std::cout << "Number_of_kd_tree_queries: "
                << number_of_point_location_queries +
                   number_of_ray_shooting_queries << std::endl;
      std::cout << "Runtime_kd_tree_queries: "
                << timer_point_location.time() +
                   timer_ray_shooting.time() << std::endl;
    }
#endif

#ifdef CGAL_NEF3_TIMER_BINARY_OPERATION
    if(cgal_nef3_timer_on) {
      timer_binary_operation.stop();
      std::cout << "Runtime_binary_operation: "
                << timer_binary_operation.time() << std::endl;
    }
#endif
    }

  };
} //namespace CGAL
#endif //CGAL_NEF_BINARY_OPERATION_H
