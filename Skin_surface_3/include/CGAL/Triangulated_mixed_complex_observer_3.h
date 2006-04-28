// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// $URL$
// $Id$
// 
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_TRIANGULATED_MIXED_COMPLEX_OBSERVER_3
#define CGAL_TRIANGULATED_MIXED_COMPLEX_OBSERVER_3

#include <CGAL/Triangulation_simplex_3.h>
#include <CGAL/Skin_surface_quadratic_surface_3.h>

CGAL_BEGIN_NAMESPACE

template <class TriangulatedMixedComplex_3,
	  class RegularTriangulation_3>
class Triangulated_mixed_complex_observer_3 {
public:
  typedef typename RegularTriangulation_3::Geom_traits     Regular_traits;
  typedef typename TriangulatedMixedComplex_3::Geom_traits Triangulated_mixed_complex_kernel;
  typedef RegularTriangulation_3                     Regular;
  typedef TriangulatedMixedComplex_3                 Triangulated_mixed_complex;

  typedef typename Regular_traits::FT                FT;

  typedef typename Regular::Vertex_handle            Rt_Vertex_handle;
  typedef typename Regular::Edge                     Rt_Edge;
  typedef typename Regular::Facet                    Rt_Facet;
  typedef typename Regular::Cell_handle              Rt_Cell_handle;
  typedef Triangulation_simplex_3<Regular>           Rt_Simplex;

  typedef typename Regular::Bare_point               Rt_Point;
  typedef typename Regular::Geom_traits              Rt_Geom_traits;
  typedef typename Rt_Geom_traits::RT                Rt_RT;
  typedef typename Regular::Weighted_point           Rt_Weighted_point;

  typedef typename Triangulated_mixed_complex::Vertex_handle TMC_Vertex_handle;
  typedef typename Triangulated_mixed_complex::Edge          TMC_Edge;
  typedef typename Triangulated_mixed_complex::Facet         TMC_Facet;
  typedef typename Triangulated_mixed_complex::Cell_handle   TMC_Cell_handle;
  
  typedef typename Triangulated_mixed_complex_kernel::Point_3        TMC_Point;
  typedef typename Triangulated_mixed_complex_kernel::RT             TMC_RT;
  typedef Weighted_point<TMC_Point,TMC_RT>        TMC_Weighted_point;
//   typedef Skin_surface_quadratic_surface_3<Polyhedron_kernel> QuadrSurface;
//   typedef Skin_surface_sphere_3<Polyhedron_kernel>         Sphere_surface;
//   typedef Skin_surface_hyperboloid_3<Polyhedron_kernel>    Hyperboloid_surface;

//   typedef typename Polyhedron_kernel::RT                   Mesh_RT;
//   typedef typename Polyhedron_kernel::Point_3              Mesh_Point;
//   typedef Weighted_point<Mesh_Point,Mesh_RT>    Mesh_Weighted_point;
  
//   typedef typename Skin_traits_3::R2P_converter R2P_converter;
//   typedef typename Skin_traits_3::T2P_converter T2P_converter;

//   Triangulated_mixed_complex_observer_3() : 
//     shrink(.5) {
//   }

  Triangulated_mixed_complex_observer_3(FT shrink) : 
    shrink(shrink) {
  }

  void after_vertex_insertion(Rt_Simplex const &sDel, 
    			      Rt_Simplex const &sVor, 
			      TMC_Vertex_handle &vh) const 
  {
  }

  void after_cell_insertion(Rt_Simplex const &s, TMC_Cell_handle &ch) const 
  {
//     if (!(s == prev_s)) {
//       prev_s = s;
//       Rt_Vertex_handle vh;
//       Rt_Edge          e;
//       Rt_Facet         f;
//       Rt_Cell_handle   ch;

//       switch (s.dimension()) {
//         case 0:
// 	  vh = s;
// 	  surf = new Sphere_surface(r2p_converter(vh->point()), shrink, 1);
// 	  break;
//         case 1:
// 	  e = s;
// 	  surf = new Hyperboloid_surface(
// 	    r2p_converter(Rt_Weighted_point(
// 	      Rt_Geom_traits().construct_weighted_circumcenter_3_object()(
// 		e.first->vertex(e.second)->point(),
// 		e.first->vertex(e.third)->point()),
// 	      Rt_Geom_traits().
// 		compute_squared_radius_smallest_orthogonal_sphere_3_object()(
// 		  e.first->vertex(e.second)->point(),
// 		  e.first->vertex(e.third)->point()))),
// 	    r2p_converter(
// 	      e.first->vertex(e.second)->point()-
// 	      e.first->vertex(e.third)->point()),
// 	    shrink, 1);
// 	  break;
//         case 2:
// 	  f = s;
// 	  surf = new Hyperboloid_surface(
// 	    r2p_converter(Rt_Weighted_point(
// 	      Rt_Geom_traits().construct_weighted_circumcenter_3_object()(
// 		f.first->vertex((f.second+1)&3)->point(),
// 		f.first->vertex((f.second+2)&3)->point(),
// 		f.first->vertex((f.second+3)&3)->point()),
// 	      Rt_Geom_traits().
// 		compute_squared_radius_smallest_orthogonal_sphere_3_object()(
// 		  f.first->vertex((f.second+1)&3)->point(),
// 		  f.first->vertex((f.second+2)&3)->point(),
// 		  f.first->vertex((f.second+3)&3)->point()))),
// 	      typename Polyhedron_kernel::Construct_orthogonal_vector_3()(
// 		r2p_converter(f.first->vertex((f.second+1)&3)->point()),
// 		r2p_converter(f.first->vertex((f.second+2)&3)->point()),
// 		r2p_converter(f.first->vertex((f.second+3)&3)->point())),
// 	    1-shrink, -1);
// 	  break;
// 	case 3:
// 	  ch = s;
// 	  surf = new Sphere_surface(
// 	    r2p_converter(Rt_Weighted_point(
// 	      Rt_Geom_traits().construct_weighted_circumcenter_3_object()(
// 		ch->vertex(0)->point(),
// 		ch->vertex(1)->point(),
// 		ch->vertex(2)->point(),
// 		ch->vertex(3)->point()),
// 	      Rt_Geom_traits().
// 		compute_squared_radius_smallest_orthogonal_sphere_3_object()(
// 		  ch->vertex(0)->point(),
// 		  ch->vertex(1)->point(),
// 		  ch->vertex(2)->point(),
// 		  ch->vertex(3)->point()))),
// 	      1-shrink, -1);
//       }

//     }
//     ch->surf = surf;
  }

  FT shrink;
  Rt_Simplex prev_s;
//   QuadrSurface *surf;
//   R2P_converter r2p_converter;
//   T2P_converter t2p_converter;
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATED_MIXED_COMPLEX_OBSERVER_3
