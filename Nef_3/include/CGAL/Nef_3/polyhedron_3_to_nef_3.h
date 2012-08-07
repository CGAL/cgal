// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// 
//
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
#ifndef CGAL_NEF_POLYHEDRON_3_TO_NEF_3_H
#define CGAL_NEF_POLYHEDRON_3_TO_NEF_3_H

#include <CGAL/Circulator_project.h>
#include <CGAL/normal_vector_newell_3.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 29
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template < class Node, class Object>
struct Project_vertex_point {
  typedef Node                  argument_type;
  typedef Object                result_type;
  Object&       operator()( Node& x)       const { return x.vertex()->point();}
  const Object& operator()( const Node& x) const { return x.vertex()->point();}
};

//SL: I added this mechanism so that it can work with
//Polyhedron_traits_with_normals_3 where the plane is
//a vector.
namespace internal
{
  template <class T>
  struct Plane_constructor;
  
  template <class K>
  struct Plane_constructor< CGAL::Plane_3<K> >
  {
    template <class Facet>
    static const CGAL::Plane_3<K>& get_plane(Facet,const CGAL::Plane_3<K>& plane){return plane;}
    static CGAL::Plane_3<K> get_type_plane(const CGAL::Point_3<K>& p,const CGAL::Vector_3<K>& vector){return CGAL::Plane_3<K>(p,vector);} 
    static CGAL::Vector_3<K> get_opposite_orthogonal_vector(const CGAL::Plane_3<K>& plane){return plane.opposite().orthogonal_vector();}
  };
  
  template <class K>
  struct Plane_constructor< CGAL::Vector_3<K> >
  {
    template <class Facet>
    static CGAL::Plane_3<K> get_plane(Facet f,const CGAL::Vector_3<K>& vector){
      return CGAL::Plane_3<K>(f->halfedge()->vertex()->point(),vector);
    }
    static const CGAL::Vector_3<K>& get_type_plane(const CGAL::Point_3<K>&,const CGAL::Vector_3<K>& vector){
      return vector;
    }
    
    static CGAL::Vector_3<K> get_opposite_orthogonal_vector(const CGAL::Vector_3<K>& vector){return -vector;}
  };

} //namespace internal

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
    CGAL_NEF_TRACEN( *point_cir);
    CGAL_NEF_TRACEN(internal::Plane_constructor<Plane>::get_type_plane(*point_cir, Vector( plane_orthogonal_vector)));
    if(plane_orthogonal_vector == Vector(0,0,0))
      std::cerr << "Error !!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    return(internal::Plane_constructor<Plane>::get_type_plane( *point_cir, Vector( plane_orthogonal_vector)));
  }
};

template<typename Items, typename Polyhedron, typename SNC_structure>
class Index_adder {
  typedef typename SNC_structure::SHalfedge_handle   SHalfedge_handle;
  
  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;
 public:
  Index_adder(Polyhedron& ) {}
  void set_hash(Halfedge_around_vertex_const_circulator,
		SHalfedge_handle) {}
  void resolve_indexes() {}
};

template<typename Polyhedron, typename SNC_structure>
class Index_adder<CGAL::SNC_indexed_items, Polyhedron, SNC_structure> {

  typedef typename SNC_structure::SHalfedge_handle   SHalfedge_handle;
  
  typedef typename Polyhedron::Halfedge_const_handle 
    Halfedge_const_handle;
  typedef typename Polyhedron::Facet_const_iterator 
    Facet_const_iterator;
  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;
  typedef typename Polyhedron::Halfedge_around_facet_const_circulator
    Halfedge_around_facet_const_circulator;
  typedef typename CGAL::Unique_hash_map<Halfedge_const_handle, 
                                         SHalfedge_handle> Hash;

  Polyhedron& P;
  Hash hash;

 public:
  Index_adder(Polyhedron& P_) : P(P_) {}

  void set_hash(Halfedge_around_vertex_const_circulator evc,
		SHalfedge_handle se) {
    hash[evc] = se;
  }
  
  void resolve_indexes() {
    Facet_const_iterator fi;
    for(fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
      Halfedge_around_facet_const_circulator 
	fc(fi->facet_begin()), end(fc);
      hash[fc]->set_index();
      hash[fc]->twin()->set_index();
      hash[fc]->twin()->source()->set_index();
      int se  = hash[fc]->get_index();
      int set = hash[fc]->twin()->get_index();
      int sv  = hash[fc]->twin()->source()->get_index();
      
      ++fc;
      CGAL_For_all(fc, end) {
	hash[fc]->set_index(se);
	hash[fc]->twin()->set_index(set);
	hash[fc]->source()->set_index(sv);
	hash[fc]->twin()->source()->set_index();
	sv = hash[fc]->twin()->source()->get_index();
      }
      hash[fc]->source()->set_index(sv);
    }
  }
};

template <class Polyhedron_, class SNC_structure>
void polyhedron_3_to_nef_3(Polyhedron_& P, SNC_structure& S)
{
  typedef Polyhedron_                                Polyhedron;
  typedef typename Polyhedron::Facet::Plane_3        Plane;
  typedef typename Polyhedron::Traits::Kernel        Kernel;
  typedef typename SNC_structure::SM_decorator       SM_decorator;
  typedef typename SNC_structure::Vertex_handle      Vertex_handle;
  typedef typename SNC_structure::SVertex_handle     SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle   SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle       SFace_handle;
  typedef typename SNC_structure::Point_3            Point_3;
  typedef typename SNC_structure::Sphere_point       Sphere_point;
  typedef typename SNC_structure::Sphere_circle      Sphere_circle;

  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator
                               Halfedge_around_vertex_const_circulator;

  Index_adder<typename SNC_structure::Items,
    Polyhedron, SNC_structure> index_adder(P);

  CGAL_NEF_TRACEN("  calculating facet's planes...");
  std::transform( P.facets_begin(), P.facets_end(),
		  P.planes_begin(), Facet_plane_3());

  typename Polyhedron::Vertex_iterator pvi;
  for( pvi = P.vertices_begin(); pvi != P.vertices_end(); ++pvi ) {
    typename Polyhedron::Vertex pv = *pvi;
    Vertex_handle nv = S.new_vertex();
    nv->point() = pv.point();
    nv->mark() = true;
    CGAL_NEF_TRACEN("v "<<pv.point());

    SM_decorator SM(&*nv);
    Halfedge_around_vertex_const_circulator pe = pv.vertex_begin(), pe_prev(pe);
    CGAL_assertion_code(Halfedge_around_vertex_const_circulator pe_0(pe));
    CGAL_assertion( pe != 0 );

    Point_3 pe_target_0(pe->opposite()->vertex()->point());
    Point_3 sp_point_0(CGAL::ORIGIN+(pe_target_0-pv.point()));
    Sphere_point sp_0(sp_point_0);
    SVertex_handle sv_0 = SM.new_svertex(sp_0);
    sv_0->mark() = true; 
    pe++;
    CGAL_assertion(pe != pv.vertex_begin());

    SVertex_handle sv_prev = sv_0;

    bool with_border = false;
    do {
      //      CGAL_assertion(!pe->is_border());
      CGAL_assertion(pe_prev->face() == pe->opposite()->face());
      CGAL_assertion(pe_prev->vertex()->point()==pv.point());
      CGAL_assertion(pe->vertex()->point()==pv.point());

      Point_3 pe_target = pe->opposite()->vertex()->point();
      Point_3 sp_point = CGAL::ORIGIN+(pe_target-pv.point());
      Sphere_point sp(sp_point);
      SVertex_handle sv = SM.new_svertex(sp);
      sv->mark() = true;
      
      //      CGAL_NEF_TRACEN(pe_prev->facet()->plane());
      CGAL_NEF_TRACEN(pe_target);
      CGAL_NEF_TRACEN(pe_prev->opposite()->vertex()->point());

      /*
      if(pe_prev->facet()->plane().is_degenerate()) {
	typename Polyhedron::Halfedge_around_facet_const_circulator fc(pv.vertex_begin()), fcend(fc);
	std::cerr << "wrong cycle "  << std::endl;
	CGAL_For_all(fc,fcend) {
	  std::cerr << "  " << fc->vertex()->point() << std::endl;
	}
      }
      */
      CGAL_assertion(pe_prev->is_border() ||
                     !internal::Plane_constructor<Plane>::get_plane(pe_prev->facet(),pe_prev->facet()->plane()).is_degenerate());
      CGAL_assertion(pe_prev->is_border() ||
		     internal::Plane_constructor<Plane>::get_plane(pe_prev->facet(),pe_prev->facet()->plane()).
		     has_on(pe_prev->opposite()->vertex()->point()));
      CGAL_assertion(pe_prev->is_border() || 
		     internal::Plane_constructor<Plane>::get_plane(pe_prev->facet(),pe_prev->facet()->plane()).has_on(pe_target));
      CGAL_assertion(pe_prev->is_border() || 
		     internal::Plane_constructor<Plane>::get_plane(pe_prev->facet(),pe_prev->facet()->plane()).has_on(pv.point()));

      if(pe_prev->is_border())
	with_border = true;
      else {
	typename Kernel::Plane_3 ss_plane
	  (CGAL::ORIGIN, 
           internal::Plane_constructor<Plane>::get_opposite_orthogonal_vector(pe_prev->facet()->plane()));
	Sphere_circle ss_circle(ss_plane);
	
	CGAL_assertion(ss_circle.has_on(sp));
	CGAL_assertion(ss_circle.has_on(sv_prev->point()));
	
	SHalfedge_handle e = SM.new_shalfedge_pair(sv_prev, sv);
	e->circle() = ss_circle;
	e->twin()->circle() = ss_circle.opposite();
	e->mark() = e->twin()->mark() = true;
	
	index_adder.set_hash(pe_prev, e);
      }

      sv_prev = sv;
      pe_prev = pe;
      ++pe;
    }
    while( pe != pv.vertex_begin() );

    CGAL_assertion(pe_prev->face() == pe_0->opposite()->face());
    CGAL_assertion(pe_prev->vertex()->point()==pv.point());
    CGAL_assertion(pe_0->vertex()->point()==pv.point());

    CGAL_NEF_TRACEN(internal::Plane_constructor<Plane>::get_plane(pe_prev->facet(),pe_prev->facet()->plane()));
    CGAL_NEF_TRACEN(pe_target_0);
    CGAL_NEF_TRACEN(pe_prev->opposite()->vertex()->point());
    CGAL_assertion(pe_prev->is_border() ||
		   !internal::Plane_constructor<Plane>::get_plane(pe_prev->facet(),pe_prev->facet()->plane()).is_degenerate());
    CGAL_assertion(pe_prev->is_border() ||
		   internal::Plane_constructor<Plane>::get_plane(pe_prev->facet(),pe_prev->facet()->plane()).
		   has_on(pe_prev->opposite()->vertex()->point()));
    CGAL_assertion(pe_prev->is_border() || 
		   internal::Plane_constructor<Plane>::get_plane(pe_prev->facet(),pe_prev->facet()->plane()).has_on(pe_target_0));
    CGAL_assertion(pe_prev->is_border() ||
		   internal::Plane_constructor<Plane>::get_plane(pe_prev->facet(),pe_prev->facet()->plane()).has_on(pv.point()));

    SHalfedge_handle e;
    if(pe_prev->is_border()) {
      with_border = true;
      e = sv_prev->out_sedge();
    } else {
      typename Kernel::Plane_3 ss_plane
	(CGAL::ORIGIN,
         internal::Plane_constructor<Plane>::get_opposite_orthogonal_vector(pe_prev->facet()->plane()));
      Sphere_circle ss_circle(ss_plane);
      
      CGAL_assertion(ss_plane.has_on(sv_prev->point()));
      CGAL_assertion(ss_circle.has_on(sp_0));
      CGAL_assertion(ss_circle.has_on(sv_prev->point()));
      
      e = SM.new_shalfedge_pair(sv_prev, sv_0);
      e->circle() = ss_circle;
      e->twin()->circle() = ss_circle.opposite();
      e->mark() = e->twin()->mark() = true;
      
      index_adder.set_hash(pe_prev, e);
    }

    // create faces
    SFace_handle fext = SM.new_sface();
    SM.link_as_face_cycle(e->twin(), fext);
    fext->mark() = false;

    if(!with_border) {
      SFace_handle fint = SM.new_sface();
      SM.link_as_face_cycle(e, fint);
      fint->mark() = false;
    }

    SM.check_integrity_and_topological_planarity();   
  }

  index_adder.resolve_indexes();
}



} //namespace CGAL

#endif //CGAL_NEF_POLYHEDRON_3_TO_NEF_3_H
