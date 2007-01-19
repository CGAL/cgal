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

#ifndef CGAL_UNION_OF_BALLS_3_H
#define CGAL_UNION_OF_BALLS_3_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// Contains the weighted converter:
#include <CGAL/Regular_triangulation_filtered_traits_3.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/triangulate_power_diagram_3.h>

// Contains the cell base of the triangulated mixed complex (and VD)
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Triangulation_simplex_3.h>

CGAL_BEGIN_NAMESPACE 

template <class UnionOfBallsTraits_3> 
class Union_of_balls_3 {
  typedef UnionOfBallsTraits_3            Gt;
  typedef Union_of_balls_3<Gt>            Self;
public:
  typedef UnionOfBallsTraits_3            Geometric_traits;
  typedef typename Gt::Weighted_point     Weighted_point;
  typedef typename Gt::Bare_point         Bare_point;
  typedef typename Gt::FT                 FT;
  // NGHK:: added for the Delaunay mesher
  typedef typename Gt::Sphere_3           Sphere;
  typedef typename Gt::Vector_3           Vector;
  
  typedef Regular_triangulation_3<Gt>     Regular;

private:
  typedef Exact_predicates_inexact_constructions_kernel     Filtered_kernel;
public:
  typedef Skin_surface_quadratic_surface_3<Filtered_kernel> 
                                                         Quadratic_surface;
public:
  typedef typename Regular::Vertex_handle                Vertex_handle;
  typedef typename Regular::Edge                         Edge;
  typedef typename Regular::Facet                        Facet;
  typedef typename Regular::Facet_circulator             Facet_circulator;
  typedef typename Regular::Cell_handle                  Cell_handle;
  typedef Triangulation_simplex_3<Regular>               Simplex;

  // pair of a del- and vor-simplex
  typedef std::pair<Simplex,Simplex>                     Anchor_point;

  //private:
  typedef typename Regular::Finite_vertices_iterator     Finite_vertices_iterator;
  typedef typename Regular::Finite_edges_iterator        Finite_edges_iterator;
  typedef typename Regular::Finite_facets_iterator       Finite_facets_iterator;
  typedef typename Regular::Finite_cells_iterator        Finite_cells_iterator;

public:
  typedef Anchor_point                                  Vertex_info;
  typedef std::pair<Simplex, Quadratic_surface *>       Cell_info;
private:
  // Triangulated_mixed_complex:
  typedef Simple_cartesian<Interval_nt_advanced>                       FK;
  typedef Triangulation_vertex_base_with_info_3<Vertex_info, FK>       Vb;
  typedef Triangulation_cell_base_with_info_3<Cell_info, FK>           Cb;
  typedef Triangulation_data_structure_3<Vb,Cb>                        Tds;
public:
  typedef Triangulation_3<FK, Tds>                                     TMC;
private:
  typedef typename TMC::Finite_vertices_iterator TMC_Vertex_iterator;
  typedef typename TMC::Finite_cells_iterator    TMC_Cell_iterator;
  typedef typename TMC::Vertex_handle            TMC_Vertex_handle;
  typedef typename TMC::Cell_handle              TMC_Cell_handle;
  typedef typename TMC::Point                    TMC_Point;
////public:
////  typedef UnionOfBallsTraits_3           Gt;
////  typedef typename Gt::Weighted_point   Weighted_point;
////  typedef typename Gt::Bare_point       Bare_point;
////  typedef typename Gt::RT               RT;
////  
////  typedef Regular_triangulation_3<Gt>     Regular;
//////   typedef Triangulation_data_structure_3 <
//////     Triangulation_vertex_base_3<GT>,
//////     Triangulated_mixed_complex_cell_3<GT, PolyhedronKernel_3> >
//////                                           Triangulated_mixed_complex_tds;
////
////  // defining the triangulated mixed complex:
////  typedef Exact_predicates_inexact_constructions_kernel    TMC_Traits;
////public:
////
////#ifdef CGAL_SKIN_SURFACE_USE_EXACT_IMPLICIT_SURFACE
////  typedef Skin_surface_quadratic_surface_3<TMC_Traits>   Quadratic_surface;
////#else
////  typedef Skin_surface_quadratic_surface_3<Simple_cartesian<double> > 
////                                                         Quadratic_surface;
////#endif // CGAL_SKIN_SURFACE_USE_EXACT_IMPLICIT_SURFACE
////
////  typedef Triangulation_3<
////    TMC_Traits,
////    Triangulation_data_structure_3
////    < Triangulated_mixed_complex_vertex_3<TMC_Traits>,
////      Triangulated_mixed_complex_cell_3<TMC_Traits,Quadratic_surface> > 
////  >                                                   Triangulated_mixed_complex;
////  typedef typename Triangulated_mixed_complex::Vertex_handle TMC_Vertex_handle;
////  typedef typename Triangulated_mixed_complex::Cell_handle   TMC_Cell_handle;
////  typedef typename TMC_Traits::Point_3                       TMC_Point;

public:  
  template < class WP_iterator >
  Union_of_balls_3(WP_iterator begin, WP_iterator end, 
		 Gt gt = Gt(),
		 bool verbose = false
		 ) 
    : gt(gt), verbose(verbose) {

    CGAL_assertion(begin != end);

    Regular regular;
    regular.insert(begin, end);
    construct_bounding_box(regular);
  
    if (verbose) {
      std::cerr << "Triangulation ready" << std::endl;
      std::cerr << "Vertices: " << regular.number_of_vertices() << std::endl;
      std::cerr << "Cells:    " << regular.number_of_cells() << std::endl;
    }
    
    // Construct the triangulated mixed complex:
    triangulate_power_diagram_3(regular, _tmc, verbose);
    
    CGAL_assertion(_tmc.is_valid());
    if (verbose) {
      std::cerr << "Triangulated mixed complex ready" << std::endl;
      std::cerr << "Vertices: " << _tmc.number_of_vertices() << std::endl;
      std::cerr << "Cells:    " << _tmc.number_of_cells() << std::endl;
    }

//     std::ofstream out("edgelength.txt");
//     for (typename Triangulated_mixed_complex::Finite_edges_iterator 
// 	   eit = _tmc.finite_edges_begin();
// 	 eit != _tmc.finite_edges_end(); eit ++) {
//       out << sqrt(_tmc.segment(eit).squared_length()) << std::endl;
//     }
  }
  const TMC &triangulated_mixed_complex() const {
    return _tmc;
  }
  
  TMC_Cell_handle locate(const TMC_Point &p) const{
    return _tmc.locate(p);
  }
private:
  void construct_bounding_box(Regular &regular);

  Gt &gt;
  TMC _tmc;
  bool verbose;
};

template <class UnionOfBallsTraits_3> 
void 
Union_of_balls_3<UnionOfBallsTraits_3>::
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
    RT dr = sqrt(CGAL::to_double(max_weight)) + 1;
  
    regular.insert(Weighted_point(
      Bare_point(bbox.xmax()+(dy+dz+dr),mid.y(),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(bbox.xmin()-(dy+dz+dr),mid.y(),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),bbox.ymax()+(dx+dz+dr),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),bbox.ymin()-(dx+dz+dr),mid.z()),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),mid.y(),bbox.zmax()+(dx+dy+dr)),-1));
    regular.insert(Weighted_point(
      Bare_point(mid.x(),mid.y(),bbox.zmin()-(dx+dy+dr)),-1));
  }
}

CGAL_END_NAMESPACE

#endif // CGAL_UNION_OF_BALLS_3_H
