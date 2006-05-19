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

#ifndef CGAL_CGAL_SKIN_SURFACE_3_H
#define CGAL_CGAL_SKIN_SURFACE_3_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// Contains the weighted converter:
#include <CGAL/Regular_triangulation_filtered_traits_3.h>

#include <CGAL/triangulate_mixed_complex_3.h>

CGAL_BEGIN_NAMESPACE 

template < class GT, 
           class QuadraticSurface, 
	   class Cb = Triangulation_cell_base_3<GT> >
class Triangulated_mixed_complex_cell_3 : public Cb
{
public:
  typedef typename Cb::Triangulation_data_structure            Triangulation_data_structure;
  typedef typename Triangulation_data_structure::Vertex_handle Vertex_handle;
  typedef typename Triangulation_data_structure::Cell_handle   Cell_handle;

  typedef QuadraticSurface                                     Quadratic_surface;
	
  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other  Cb2;
    typedef Triangulated_mixed_complex_cell_3<GT, Quadratic_surface, Cb2>
                                                           Other;
  };

  Triangulated_mixed_complex_cell_3() : Cb() {
  }
  Triangulated_mixed_complex_cell_3(Vertex_handle v0, Vertex_handle v1,
				    Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {
  }

  Quadratic_surface *surf;
};

template <class SkinSurfaceTraits_3> 
class Skin_surface_3 {
public:
  typedef SkinSurfaceTraits_3           Gt;
  typedef typename Gt::Weighted_point   Weighted_point;
  typedef typename Gt::Bare_point       Bare_point;
  typedef typename Gt::FT               FT;
  
  typedef Regular_triangulation_3<Gt>     Regular;
//   typedef Triangulation_data_structure_3 <
//     Triangulation_vertex_base_3<GT>,
//     Triangulated_mixed_complex_cell_3<GT, PolyhedronKernel_3> >
//                                           Triangulated_mixed_complex_tds;

  // defining the triangulated mixed complex:
  typedef Exact_predicates_inexact_constructions_kernel    TMC_Traits;
public:

#ifdef CGAL_SKIN_SURFACE_USE_EXACT_IMPLICIT_SURFACE
  typedef Skin_surface_quadratic_surface_3<TMC_Traits>   Quadratic_surface;
#else
  typedef Skin_surface_quadratic_surface_3<Simple_cartesian<double> > 
                                                         Quadratic_surface;
#endif // CGAL_SKIN_SURFACE_USE_EXACT_IMPLICIT_SURFACE

  typedef Triangulation_3<
    TMC_Traits,
    Triangulation_data_structure_3
    < Triangulation_vertex_base_3<TMC_Traits>,
      Triangulated_mixed_complex_cell_3<TMC_Traits,Quadratic_surface> > 
  >                                                   Triangulated_mixed_complex;
  typedef typename Triangulated_mixed_complex::Vertex_handle TMC_Vertex_handle;
  typedef typename Triangulated_mixed_complex::Cell_handle   TMC_Cell_handle;
  typedef typename TMC_Traits::Point_3                       TMC_Point;
  
  template < class WP_iterator >
  Skin_surface_3(WP_iterator begin, WP_iterator end, 
		 FT shrink_factor,
		 bool grow_balls = true,
		 Gt gt = Gt(),
		 bool verbose = false
		 ) 
    : gt(gt), shrink(shrink_factor), verbose(verbose) {

    CGAL_assertion(begin != end);

    Regular regular;
    if (grow_balls) {
      for (; begin != end; begin++) {
	regular.insert(Weighted_point(*begin, begin->weight()/shrink));
      }
    } else {
      regular.insert(begin, end);
    }
    construct_bounding_box(regular);
  
    if (verbose) {
      std::cerr << "Triangulation ready" << std::endl;
      std::cerr << "Vertices: " << regular.number_of_vertices() << std::endl;
      std::cerr << "Cells:    " << regular.number_of_cells() << std::endl;
    }
    
    // Construct the triangulated mixed complex:
    triangulate_mixed_complex_3(regular, shrink, _tmc,verbose);
    
    CGAL_assertion(_tmc.is_valid());
    if (verbose) {
      std::cerr << "Triangulated mixed complex ready" << std::endl;
      std::cerr << "Vertices: " << _tmc.number_of_vertices() << std::endl;
      std::cerr << "Cells:    " << _tmc.number_of_cells() << std::endl;
    }
//     std::ofstream out("vertices.txt");
//     for (typename Triangulated_mixed_complex::Finite_vertices_iterator 
// 	   vit = _tmc.finite_vertices_begin();
// 	 vit != _tmc.finite_vertices_end(); vit ++) {
//       out << vit->point().x().exact() << std::endl;
//     }
  }
  const Triangulated_mixed_complex &triangulated_mixed_complex() const {
    return _tmc;
  }
  
  TMC_Cell_handle locate(const TMC_Point &p) const{
    return _tmc.locate(p);
  }
private:
  void construct_bounding_box(Regular &regular);

  Gt &gt;
  FT shrink;
  Triangulated_mixed_complex _tmc;
  bool verbose;
};

template <class SkinSurfaceTraits_3> 
void 
Skin_surface_3<SkinSurfaceTraits_3>::
construct_bounding_box(Regular &regular) 
{
  typedef typename Regular::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Regular::Geom_traits     GT;
  typedef typename GT::Bare_point             Point;
  typedef typename GT::Point                Weighted_point;
  typedef typename GT::RT                     RT;
  
  Finite_vertices_iterator vit = regular.finite_vertices_begin();
  if (vit != regular.finite_vertices_end()) {
    Bbox_3 bbox = vit->point().bbox();
    RT max_weight=vit->point().weight();
    while (++vit != regular.finite_vertices_end()) {
      bbox = bbox + vit->point().bbox();
      if (max_weight < vit->point().weight())
	max_weight = vit->point().weight();
    }

    // add a bounding octahedron:
    RT dx = bbox.xmax() - bbox.xmin();
    RT dy = bbox.ymax() - bbox.ymin();
    RT dz = bbox.zmax() - bbox.zmin();
  
    Bare_point mid(bbox.xmin() + dx/2, bbox.ymin() + dy/2, bbox.zmin() + dz/2);
    RT dr = sqrt(CGAL::to_double(max_weight)) + .001;
  
    regular.insert(Weighted_point(
      Bare_point(bbox.xmax()+(dy+dz+dr)/shrink,mid.y(),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(bbox.xmin()-(dy+dz+dr)/shrink,mid.y(),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),bbox.ymax()+(dx+dz+dr)/shrink,mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),bbox.ymin()-(dx+dz+dr)/shrink,mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),mid.y(),bbox.zmax()+(dx+dy+dr)/shrink),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),mid.y(),bbox.zmin()-(dx+dy+dr)/shrink),-1));
  }
}

// template <class InputIterator, class Polyhedron_3, class SkinSurfaceTraits_3>
// void skin_surface_3(InputIterator first, InputIterator last,
//   Polyhedron_3 &polyhedron, const SkinSurfaceTraits_3 &skin_surface_traits,
//   bool verbose = false) {
//   if (first == last) {
//     return;
//   }

//   // Types
//   typedef SkinSurfaceTraits_3                              Skin_surface_traits;
//   typedef typename Skin_surface_traits::Regular_traits     Regular_traits;
//   typedef typename Regular_traits::Bare_point              Reg_point;
//   typedef typename Regular_traits::Weighted_point          Reg_weighted_point;

//   typedef Regular_triangulation_3<Regular_traits> Regular;
//   typedef Triangulated_mixed_complex_3<SkinSurfaceTraits_3>
//                                                   Triangulated_mixed_complex;
//   typedef Marching_tetrahedra_traits_skin_surface_3<
//     Triangulated_mixed_complex,
//     Polyhedron_3,
//     typename SkinSurfaceTraits_3::T2P_converter>  Marching_tetrahedra_traits;
//   typedef Marching_tetrahedra_observer_default_3<
//     Triangulated_mixed_complex, Polyhedron_3>     Marching_tetrahedra_observer;
    
//   // Code
//   Regular regular;
//   Triangulated_mixed_complex triangulated_mixed_complex;

//   while (first != last) {
//     regular.insert((*first));
//     first++;
//   }

//   skin_surface_construct_bounding_box_3(regular,skin_surface_traits);
  
//   if (verbose) {
//     std::cerr << "Triangulation ready" << std::endl;
//   }

//   // Construct the triangulated mixed complex:
//   triangulate_mixed_complex_3(
//     regular, triangulated_mixed_complex, skin_surface_traits);

//   CGAL_assertion(triangulated_mixed_complex.is_valid());
//   if (verbose) {
//     std::cerr << "Triangulated mixed complex ready" << std::endl;
//   }

//   // Extract the coarse mesh using marching_tetrahedra
//   Marching_tetrahedra_traits marching_traits;
//   marching_tetrahedra_3(
//     triangulated_mixed_complex, polyhedron, marching_traits);

//   if (verbose) {
//     std::cerr << "Mesh ready" << std::endl;
//   }
  
// }

CGAL_END_NAMESPACE

#endif // CGAL_CGAL_SKIN_SURFACE_TRAITS_3_H
