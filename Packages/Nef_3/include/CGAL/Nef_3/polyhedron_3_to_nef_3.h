//#define CGAL_P2NEF3_USE_SM_OVERLAY
// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/polyhedron_3_to_nef_3.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// Constructor of a SNC structure from a Polyhedron
// ============================================================================
#ifndef POLYHEDRON_3_TO_NEF_3_H
#define POLYHEDRON_3_TO_NEF_3_H

#include <CGAL/Circulator_project.h>
#include <CGAL/normal_vector_newell_3.h>

#undef _DEBUG
#define _DEBUG 5
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template < class Node, class Object>
struct Project_vertex_point {
  typedef Node                  argument_type;
  typedef Object                result_type;
  Object&       operator()( Node& x)       const { return x.vertex()->point();}
  const Object& operator()( const Node& x) const { return x.vertex()->point();}
};

struct Facet_plane_3 {
  template < class Facet_>
  typename Facet_::Plane_3 operator()( Facet_& f) {
    typedef Facet_                              Facet;
    typedef typename Facet::Plane_3             Plane;
    typedef Kernel_traits< Plane>               KernelTraits;
    typedef typename KernelTraits::Kernel       Kernel;
    typedef typename Kernel::Vector_3           Vector;
    typedef typename Kernel::Direction_3        Direction;
    typedef typename Facet::Halfedge_around_facet_const_circulator
                                                Halfedge_circulator;
    typedef typename Facet::Halfedge            Halfedge;
    typedef typename Halfedge::Vertex           Vertex;
    typedef typename Vertex::Point_3            Point;
    typedef Project_vertex_point< Halfedge, const Point> Proj_vertex_point;
    typedef Circulator_project< Halfedge_circulator, Proj_vertex_point,
      const Point, const Point*> Circulator;
    /* TODO: to implement a better approach
       typedef Project_vertex< Halfedge> Project_vertex;
       typedef Project_point< Vertex> Project_point;
       typedef Compose< Project_vertex, Project_point> Projector;
       typedef Circulator_project< Halfedge_circulator, Projector> Circulator;
    */
    Circulator point_cir( f.facet_begin());
    Vector plane_orthogonal_vector;
    normal_vector_newell_3( point_cir, point_cir, plane_orthogonal_vector);
    return( Plane( *point_cir, Direction( plane_orthogonal_vector)));
  }
};


template <class Polyhedron_, class SNC_structure, class SNC_constructor>
void polyhedron_3_to_nef_3(Polyhedron_& P, SNC_structure& S)
{
  typedef Polyhedron_                                Polyhedron;
  typedef typename SNC_constructor::SM_overlayer     SM_overlayer;
  typedef typename SNC_structure::SNC_decorator      SNC_decorator;
  typedef typename SNC_structure::SM_decorator       SM_decorator;
  typedef typename SNC_structure::Vertex_handle      Vertex_handle;
  typedef typename SNC_structure::SVertex_handle     SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle   SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle       SFace_handle;
  typedef typename SNC_structure::Point_3            Point_3;
  typedef typename SNC_structure::Sphere_point       Sphere_point;
  typedef typename SNC_structure::Sphere_segment     Sphere_segment;
  typedef typename SNC_structure::Sphere_circle      Sphere_circle;
                                  
  std::transform( P.facets_begin(), P.facets_end(),
		  P.planes_begin(), Facet_plane_3());
  /* determine the plane of polyhedron's facet */
  // CGAL::set_pretty_mode( std::cout);
  // std::copy( P.planes_begin(), P.planes_end(), std::ostream_iterator
  //     < typename Polyhedron::Plane_3>( std::cout, "\n"));

  SNC_decorator D(S);
  SNC_constructor C(S);
    
  typename Polyhedron::Vertex_iterator pvi;
  for( pvi = P.vertices_begin(); pvi != P.vertices_end(); ++pvi ) {
    typename Polyhedron::Vertex pv = *pvi;
    Vertex_handle nv = S.new_vertex();
    D.point(nv) = pv.point();
    D.mark(nv) = true;
    TRACEN("v "<<pv.point());

#ifdef CGAL_P2NEF3_USE_SM_OVERLAY
    std::list<Sphere_point> sp_list;
    std::list<Sphere_segment> ss_list;
    typename Polyhedron::Halfedge_around_vertex_const_circulator 
      ph_edge = pv.vertex_begin();
    if(ph_edge != 0) {
	do {
	  Point_3 
	    ph_point(ph_edge++->opposite()->vertex()->point());
	  Sphere_point sp(ph_point-pv.point());
	  sp_list.push_back(sp);
	}
	while(ph_edge != pv.vertex_begin());
	
	std::list<Sphere_point>::iterator sp_prev, sp_cur;
	ph_edge = pv.vertex_begin();
	sp_cur = sp_list.begin();
	sp_prev = sp_cur++;
	while(sp_cur != sp_list.end()) {
	    Sphere_circle
	      ss_circle(ph_edge++->facet()->plane().opposite());
	    CGAL_assertion(ss_circle.has_on(*sp_prev));
	    CGAL_assertion(ss_circle.has_on(*sp_cur));
	    ss_list.
	      push_back(Sphere_segment(*sp_prev++,*sp_cur++,ss_circle));
	  }
	Sphere_circle
	  ss_circle(ph_edge++->facet()->plane().opposite());
	CGAL_assertion(ss_circle.has_on(*sp_prev));
	CGAL_assertion(ss_circle.has_on(*(sp_list.begin())));
	ss_list.
	  push_back(Sphere_segment(*sp_prev,*(sp_list.begin()),ss_circle));
	CGAL_assertion(ph_edge == pv.vertex_begin());
      }
    SM_overlayer O(nv);
    O.create_from_segments(ss_list.begin(),ss_list.end());
    O.simplify();
    
#else // CGAL_P2NEF3_USE_SM_OVERLAY

    SM_decorator SM(nv);
    typename Polyhedron::Halfedge_around_vertex_const_circulator
      pe = pv.vertex_begin(), pe_prev(pe), pe_0(pe);
    CGAL_assertion( pe != 0 );

    Point_3 pe_target_0(pe->opposite()->vertex()->point());
    Point_3 sp_point_0(CGAL::ORIGIN+(pe_target_0-pv.point()));
    Sphere_point sp_0(sp_point_0);
    SVertex_handle sv_0 = SM.new_vertex(sp_0);
    SM.mark(sv_0) = true; 
    pe++;
    CGAL_assertion(pe != pv.vertex_begin());

    SVertex_handle sv_prev = sv_0;

    do {
      CGAL_assertion(!pe->is_border());
      CGAL_assertion(pe_prev->face() == pe->opposite()->face());
      CGAL_assertion(pe_prev->vertex()->point()==pv.point());
      CGAL_assertion(pe->vertex()->point()==pv.point());

      Point_3 pe_target = pe->opposite()->vertex()->point();
      Point_3 sp_point = CGAL::ORIGIN+(pe_target-pv.point());
      Sphere_point sp(sp_point);
      SVertex_handle sv = SM.new_vertex(sp);
      SM.mark(sv) = true;
      
      TRACEN(pe_prev->facet()->plane());
      TRACEN(pe_target);
      TRACEN(pe_prev->opposite()->vertex()->point());
      CGAL_assertion(!pe_prev->facet()->plane().is_degenerate());
      CGAL_assertion(pe_prev->facet()->plane().
		     has_on(pe_prev->opposite()->vertex()->point()));
      CGAL_assertion(pe_prev->facet()->plane().has_on(pe_target));
      CGAL_assertion(pe_prev->facet()->plane().has_on(pv.point()));

      typename Polyhedron::Facet::Plane_3 ss_plane
	(CGAL::ORIGIN, 
	 pe_prev->facet()->plane().opposite().orthogonal_vector());
      Sphere_circle ss_circle(ss_plane);

      //CGAL_assertion(!ss_plane.is_degenerate());
      //CGAL_assertion(ss_plane.has_on(CGAL::ORIGIN));
      //CGAL_assertion(ss_plane.has_on(sp_point));
      //CGAL_assertion(ss_plane.has_on(SM.point(sv_prev));
      CGAL_assertion(ss_circle.has_on(sp));
      CGAL_assertion(ss_circle.has_on(SM.point(sv_prev)));

      SHalfedge_handle e = SM.new_edge_pair(sv_prev, sv);
      SM.circle(e) = ss_circle;
      SM.circle(SM.twin(e)) = ss_circle.opposite();
      SM.mark(e) = true;
	  
      sv_prev = sv;
      pe_prev = pe;
      pe++;
    }
    while( pe != pv.vertex_begin() );

    CGAL_assertion(!pe_0->is_border());
    CGAL_assertion(pe_prev->face() == pe_0->opposite()->face());
    CGAL_assertion(pe_prev->vertex()->point()==pv.point());
    CGAL_assertion(pe_0->vertex()->point()==pv.point());
    
    TRACEN(pe_prev->facet()->plane());
    TRACEN(pe_target_0);
    TRACEN(pe_prev->opposite()->vertex()->point());
    CGAL_assertion(!pe_prev->facet()->plane().is_degenerate());
    CGAL_assertion(pe_prev->facet()->plane().
		   has_on(pe_prev->opposite()->vertex()->point()));
    CGAL_assertion(pe_prev->facet()->plane().has_on(pe_target_0));
    CGAL_assertion(pe_prev->facet()->plane().has_on(pv.point()));

    typename Polyhedron::Facet::Plane_3 ss_plane
      (CGAL::ORIGIN,
       pe_prev->facet()->plane().opposite().orthogonal_vector());
    Sphere_circle ss_circle(ss_plane);

    //CGAL_assertion(!ss_plane.is_degenerate());
    //CGAL_assertion(ss_plane.has_on(CGAL::ORIGIN));
    //CGAL_assertion(ss_plane.has_on(sp_point_0));
    CGAL_assertion(ss_plane.has_on(SM.point(sv_prev)));
    CGAL_assertion(ss_circle.has_on(sp_0));
    CGAL_assertion(ss_circle.has_on(SM.point(sv_prev)));
    
    SHalfedge_handle e = SM.new_edge_pair(sv_prev, sv_0);
    SM.circle(e) = ss_circle;
    SM.circle(SM.twin(e)) = ss_circle.opposite();
    SM.mark(e) = true;

    // create facets
    SFace_handle fint = SM.new_face();
    SFace_handle fext = SM.new_face();
    SM.link_as_face_cycle(e, fint);
    SM.link_as_face_cycle(SM.twin(e), fext);
    
    // mark facets properly...
    SM.mark(fint) = true;
    SM.mark(fext) = false;
    SM.check_integrity_and_topological_planarity();
    
    /* Why not? SM_overlayer O(nv); O.simplify(); */

#endif // CGAL_P2NEF3_USE_SM_OVERLAY
	
#ifdef SM_VISUALIZE
    CGAL::OGL::add_sphere();
    SMV V(nv, CGAL::OGL::spheres_.back());
    V.draw_map();
#endif
  }

#ifdef SM_VISUALIZE
  CGAL::OGL::start_viewer();
#endif

  return;
}

CGAL_END_NAMESPACE

#endif //POLYHEDRON_3_TO_NEF_3_H

