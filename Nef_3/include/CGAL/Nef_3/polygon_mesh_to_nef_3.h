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
//                 Andreas Fabri

#ifndef CGAL_POLYGON_MESH_TO_NEF_3_H
#define CGAL_POLYGON_MESH_TO_NEF_3_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/normal_vector_newell_3.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/boost/graph/helpers.h>

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

template<typename Items, typename Polyhedron, typename SNC_structure, typename HalfedgeIndexMap>
class Face_graph_index_adder {
  typedef typename SNC_structure::SHalfedge_handle   SHalfedge_handle;

  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
 public:
  Face_graph_index_adder(Polyhedron&, HalfedgeIndexMap ) {}
  void set_hash(halfedge_descriptor,
		SHalfedge_handle) {}
  void resolve_indexes() {}
};

template<typename PolygonMesh, typename SNC_structure, typename HalfedgeIndexMap>
class Face_graph_index_adder<CGAL::SNC_indexed_items, PolygonMesh, SNC_structure, HalfedgeIndexMap> {

  typedef typename SNC_structure::SHalfedge_handle   SHalfedge_handle;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

  typedef Halfedge_around_face_circulator<PolygonMesh>
    Halfedge_around_facet_const_circulator;
    typedef std::vector<SHalfedge_handle> Hash;

  PolygonMesh& P;
  HalfedgeIndexMap him;
  Hash hash;

public:
  Face_graph_index_adder(PolygonMesh& P_, HalfedgeIndexMap him) : P(P_), him(him)
  {
    hash.resize(num_halfedges(P));
  }

  void set_hash(halfedge_descriptor evc,
		SHalfedge_handle se) {
    hash[get(him,evc)] = se;
  }
  
  void resolve_indexes()
  {
    BOOST_FOREACH(face_descriptor fi, faces(P)) {
      Halfedge_around_facet_const_circulator 
	fc(halfedge(fi,P),P), end(fc);
      typename boost::property_traits<HalfedgeIndexMap>::value_type
        index = get(him,*fc);
      hash[index]->set_index();
      hash[index]->twin()->set_index();
      hash[index]->twin()->source()->set_index();
      int se  = hash[index]->get_index();
      int set = hash[index]->twin()->get_index();
      int sv  = hash[index]->twin()->source()->get_index();
      
      ++fc;
      CGAL_For_all(fc, end) {
        index = get(him,*fc);
	hash[index]->set_index(se);
	hash[index]->twin()->set_index(set);
	hash[index]->source()->set_index(sv);
	hash[index]->twin()->source()->set_index();
	sv = hash[index]->twin()->source()->get_index();
      }
      hash[get(him,*fc)]->source()->set_index(sv);
    }
  }
};

template <class PolygonMesh, class SNC_structure, class FaceIndexMap, class HalfedgeIndexMap>
void polygon_mesh_to_nef_3(PolygonMesh& P, SNC_structure& S, FaceIndexMap fimap, HalfedgeIndexMap himap)
{
  typedef typename boost::property_map<PolygonMesh, vertex_point_t>::type PMap;
  typedef typename SNC_structure::Plane_3                   Plane;
  typedef typename SNC_structure::Vector_3                           Vector_3;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename SNC_structure::SM_decorator       SM_decorator;
  typedef typename SNC_structure::Vertex_handle      Vertex_handle;
  typedef typename SNC_structure::SVertex_handle     SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle   SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle       SFace_handle;
  typedef typename SNC_structure::Point_3            Point_3;
  typedef typename SNC_structure::Sphere_point       Sphere_point;
  typedef typename SNC_structure::Sphere_circle      Sphere_circle;

  typedef Halfedge_around_target_circulator<PolygonMesh>
                               Halfedge_around_vertex_const_circulator;


  PMap pmap = get(CGAL::vertex_point,P);

  std::vector<Vector_3> normals(num_faces(P));

  BOOST_FOREACH(face_descriptor f, faces(P)){
    Vertex_around_face_circulator<PolygonMesh> vafc(halfedge(f,P),P), done(vafc);
    Vector_3 v;
    normal_vector_newell_3(vafc, done, pmap, v);
    normals[get(fimap,f)] = - v;
  }

  Face_graph_index_adder<typename SNC_structure::Items,
                 PolygonMesh, SNC_structure,HalfedgeIndexMap> index_adder(P,himap);


  BOOST_FOREACH(vertex_descriptor pv, vertices(P) ) {
    Vertex_handle nv = S.new_vertex();
    nv->point() = get(pmap,pv);
    nv->mark() = true;
    CGAL_NEF_TRACEN("v "<< get(pmap,pv));

    SM_decorator SM(&*nv);
    Halfedge_around_vertex_const_circulator pec(pv,P), pec_prev(pec), done(pec);
    halfedge_descriptor pe = *pec, pe_prev= *pec_prev;
    CGAL_assertion_code(Halfedge_around_vertex_const_circulator pe_0(pec));
    // CGAL_assertion( pe != 0 );

    Point_3 pe_target_0(get(pmap,target(opposite(pe,P),P)));
    Point_3 sp_point_0(CGAL::ORIGIN+(pe_target_0 - get(pmap,pv)));
    Sphere_point sp_0(sp_point_0);
    SVertex_handle sv_0 = SM.new_svertex(sp_0);
    sv_0->mark() = true; 
    pec++;
    pe = *pec;
    //CGAL_assertion(pe != pv->vertex_begin());

    SVertex_handle sv_prev = sv_0;

    bool with_border = false;
    do {
      CGAL_assertion(face(pe_prev,P) == face(opposite(pe,P),P));
      CGAL_assertion(get(pmap,target(pe_prev,P)) == get(pmap,pv));
      CGAL_assertion(get(pmap,target(pe,P)) == get(pmap,pv));

      Point_3 pe_target = get(pmap,target(opposite(pe,P),P));
      Point_3 sp_point = CGAL::ORIGIN+(pe_target - get(pmap,pv));
      Sphere_point sp(sp_point);
      SVertex_handle sv = SM.new_svertex(sp);
      sv->mark() = true;
      
      CGAL_NEF_TRACEN(pe_target);
      CGAL_NEF_TRACEN(get(pmap,target(opposite(pe_prev,P),P)));

      if(is_border(pe_prev,P))
	with_border = true;
      else {
  Plane ss_plane( CGAL::ORIGIN, normals[get(fimap,face(pe_prev,P))] );
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
      pec_prev = pec;
      ++pec;
      pe = *pec;
      pe_prev = *pec_prev;
    }
    while( pec != done );

    CGAL_assertion(face(pe_prev,P) == face(opposite(*pe_0,P),P));
    CGAL_assertion(get(pmap,target(pe_prev,P)) == get(pmap,pv));
    CGAL_assertion(get(pmap,target(*pe_0,P)) == get(pmap,pv));

    SHalfedge_handle e;
    if(is_border(pe_prev,P)) {
      with_border = true;
      e = sv_prev->out_sedge();
    } else {
      Plane ss_plane( CGAL::ORIGIN, normals[get(fimap,face(pe_prev,P))] );
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

template <class Polyhedron, class SNC_structure>
void polyhedron_3_to_nef_3(Polyhedron& P, SNC_structure& S)
{
  typedef typename boost::property_map<Polyhedron, face_external_index_t>::type FIMap;
  FIMap fimap = get(CGAL::face_external_index,P);
  typedef typename boost::property_map<Polyhedron, halfedge_external_index_t>::type HIMap;
  HIMap himap = get(CGAL::halfedge_external_index,P);
  polygon_mesh_to_nef_3(P, S, fimap, himap);
}

template <class SM, class SNC_structure>
void polygon_mesh_to_nef_3(SM& sm, SNC_structure& snc)
{
  typedef typename boost::property_map<SM, face_index_t>::type FIMap;
  FIMap fimap = get(CGAL::face_index,sm);
  typedef typename boost::property_map<SM, boost::halfedge_index_t>::type HIMap;
  HIMap himap = get(boost::halfedge_index,sm);

  polygon_mesh_to_nef_3(sm, snc, fimap, himap);
}


} //namespace CGAL

#endif //CGAL_POLYGON_MESH_TO_NEF_3_H
