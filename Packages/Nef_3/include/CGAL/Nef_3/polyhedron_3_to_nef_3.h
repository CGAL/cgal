// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
#ifndef POLYHEDRON_3_TO_NEF_3_H
#define POLYHEDRON_3_TO_NEF_3_H

#include <CGAL/Circulator_project.h>
#include <CGAL/normal_vector_newell_3.h>
#include <CGAL/Nef_S2/SM_point_locator.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 29
#include <CGAL/Nef_2/debug.h>

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
    TRACEN( *point_cir);
    TRACEN(Plane( *point_cir, Vector( plane_orthogonal_vector)));
    return( Plane( *point_cir, Vector( plane_orthogonal_vector)));
  }
};


template <class Polyhedron_, class SNC_structure>
void polyhedron_3_to_nef_3(Polyhedron_& P, SNC_structure& S)
{
  typedef Polyhedron_                                Polyhedron;
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

  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator
                               Halfedge_around_vertex_const_circulator;
                                  
  TRACEN("  calculating facet's planes...");
  std::transform( P.facets_begin(), P.facets_end(),
		  P.planes_begin(), Facet_plane_3());
  SNC_decorator D(S);

  //  Progress_indicator_clog progress
  //    ( P.size_of_vertices(),
  //      "polyhedron_3_to_nef_3: constructing local view of vertices...");

  typename Polyhedron::Vertex_iterator pvi;
  for( pvi = P.vertices_begin(); pvi != P.vertices_end(); ++pvi ) {
    //    progress++;
    typename Polyhedron::Vertex pv = *pvi;
    Vertex_handle nv = S.new_vertex();
    D.point(nv) = pv.point();
    D.mark(nv) = true;
    TRACEN("v "<<pv.point());

    SM_decorator SM(&*nv);
    Halfedge_around_vertex_const_circulator pe = pv.vertex_begin(), pe_prev(pe);
    CGAL_assertion_code(Halfedge_around_vertex_const_circulator pe_0(pe));
    CGAL_assertion( pe != 0 );

    Point_3 pe_target_0(pe->opposite()->vertex()->point());
    Point_3 sp_point_0(CGAL::ORIGIN+(pe_target_0-pv.point()));
    Sphere_point sp_0(sp_point_0);
    SVertex_handle sv_0 = SM.new_svertex(sp_0);
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
      SVertex_handle sv = SM.new_svertex(sp);
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

      CGAL_assertion(ss_circle.has_on(sp));
      CGAL_assertion(ss_circle.has_on(SM.point(sv_prev)));

      SHalfedge_handle e = SM.new_shalfedge_pair(sv_prev, sv);
      SM.circle(e) = ss_circle;
      SM.circle(SM.twin(e)) = ss_circle.opposite();
      SM.mark(e) = SM.mark(SM.twin(e)) = true;
	  
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

    CGAL_assertion(ss_plane.has_on(SM.point(sv_prev)));
    CGAL_assertion(ss_circle.has_on(sp_0));
    CGAL_assertion(ss_circle.has_on(SM.point(sv_prev)));
    
    SHalfedge_handle e = SM.new_shalfedge_pair(sv_prev, sv_0);
    SM.circle(e) = ss_circle;
    SM.circle(SM.twin(e)) = ss_circle.opposite();
    SM.mark(e) = SM.mark(SM.twin(e)) = true;

    // create faces
    SFace_handle fint = SM.new_sface();
    SFace_handle fext = SM.new_sface();
    SM.link_as_face_cycle(e, fint);
    SM.link_as_face_cycle(SM.twin(e), fext);
    
    // mark faces properly...
    SM.mark(fint) = true;
    SM.mark(fext) = false;
    SM.check_integrity_and_topological_planarity();   

    //    SM_point_locator L(nv);
    //    L.init_marks_of_halfspheres(); 
  }
}

CGAL_END_NAMESPACE

#endif //POLYHEDRON_3_TO_NEF_3_H

