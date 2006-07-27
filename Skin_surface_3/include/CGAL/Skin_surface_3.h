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

#ifndef CGAL_SKIN_SURFACE_3_H
#define CGAL_SKIN_SURFACE_3_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// Contains the weighted converter:
#include <CGAL/Regular_triangulation_filtered_traits_3.h>

#include <CGAL/triangulate_mixed_complex_3.h>

// Needed for the (Delaunay) surface mesher
#include <CGAL/Skin_surface_mesher_oracle_3.h>
#include <CGAL/Triangulation_simplex_3.h>

CGAL_BEGIN_NAMESPACE 

template < class GT, 
           class SkinSurface_3, 
	   class Cb = Triangulation_cell_base_3<GT> >
class Triangulated_mixed_complex_cell_3 : public Cb
{
public:
  typedef typename Cb::Triangulation_data_structure            Triangulation_data_structure;
  typedef typename Triangulation_data_structure::Vertex_handle Vertex_handle;
  typedef typename Triangulation_data_structure::Cell_handle   Cell_handle;

  typedef typename SkinSurface_3::Quadratic_surface           Quadratic_surface;
  typedef typename SkinSurface_3::Simplex                     Simplex;
	
  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other  Cb2;
    typedef Triangulated_mixed_complex_cell_3<GT, SkinSurface_3, Cb2>
                                                           Other;
  };

  Triangulated_mixed_complex_cell_3() : Cb() {
  }
  Triangulated_mixed_complex_cell_3(Vertex_handle v0, Vertex_handle v1,
				    Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {
  }

//   template <class Input_point>
//   Sign sign(const Input_point &p) const {
//     return surf->sign(p);
//   }
  Quadratic_surface *surf;
  Simplex simp;
};

template < class GT, 
	   class Vb = Triangulation_vertex_base_3<GT> >
class Triangulated_mixed_complex_vertex_3 : public Vb
{
public:
  typedef typename Vb::Point           Point;
  typedef typename Vb::Cell_handle     Cell_handle;

  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef Triangulated_mixed_complex_vertex_3<GT, Vb2>   Other;
  };

  Triangulated_mixed_complex_vertex_3() {}
  Triangulated_mixed_complex_vertex_3(const Point&p)                : Vb(p) {}
  Triangulated_mixed_complex_vertex_3(const Point&p, Cell_handle c) : Vb(p, c) {}

  Sign sign() const {
    return Vb::cell()->surf->sign(Vb::point());
  }
};

template <class SkinSurfaceTraits_3> 
class Skin_surface_3 {
  typedef SkinSurfaceTraits_3             Gt;
  typedef Skin_surface_3<Gt>              Self;
public:
  typedef SkinSurfaceTraits_3             Geometric_traits;
  typedef typename Gt::Weighted_point     Weighted_point;
  typedef typename Weighted_point::Weight RT;
  // NGHK:: added for the Delaunay mesher
  typedef typename Gt::Sphere_3           Sphere_3;
private:
  typedef typename Weighted_point::Point  Bare_point;
  
  typedef Regular_triangulation_3<Gt>     Regular;

public:
  // NGHK: remove later?
  typedef Triangulation_simplex_3<Regular>               Simplex;

  // defining the triangulated mixed complex:
  typedef Exact_predicates_exact_constructions_kernel    TMC_traits;

  typedef Skin_surface_quadratic_surface_3<TMC_traits>   Quadratic_surface;

  typedef Triangulation_3<
    TMC_traits,
    Triangulation_data_structure_3
    < Triangulated_mixed_complex_vertex_3<TMC_traits>,
      Triangulated_mixed_complex_cell_3<TMC_traits,Self> > 
  >                                      Triangulated_mixed_complex;

  typedef typename Triangulated_mixed_complex::Vertex_handle TMC_Vertex_handle;
  typedef typename Triangulated_mixed_complex::Cell_handle   TMC_Cell_handle;

  // NGHK: added for the (Delaunay) surface mesher, document
  typedef Exact_predicates_inexact_constructions_kernel  Mesher_Gt;
  typedef Skin_surface_mesher_oracle_3<Mesher_Gt,Self> Surface_mesher_traits_3;

private:
  typedef typename TMC_traits::Point_3                       TMC_Point;
  
public:
  template < class WP_iterator >
  Skin_surface_3(WP_iterator begin, WP_iterator end, 
		 RT shrink_factor,
		 bool grow_balls = true,
		 Gt gt = Gt(),
		 bool _verbose = false
		 ) 
    : gt(gt), shrink(shrink_factor), verbose(_verbose) {

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
    triangulate_mixed_complex_3(regular, shrink, _tmc, verbose);
    
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
  const Triangulated_mixed_complex &triangulated_mixed_complex() const {
    return _tmc;
  }
  
  TMC_Cell_handle locate(const TMC_Point &p) const{
    last_ch = _tmc.locate(p, last_ch);
    return last_ch;
  }


  // NGHK: added for the (Delaunay) surface mesher, document
  Sphere_3 bounding_sphere() const {
    return _bounding_sphere;
  }
  RT squared_error_bound() const {
    return .01;
  }
  Sign operator()(const Bare_point &p) const {
     Cartesian_converter<typename Bare_point::R, TMC_traits > converter;
     TMC_Point p_tmc = converter(p);
     TMC_Cell_handle ch = locate(p_tmc);
     if (_tmc.is_infinite(ch)) {
       // Infinite cells do not have a pointer to a surface
       return NEGATIVE;
     }
     return ch->surf->sign(p_tmc);
  }
  typename Mesher_Gt::FT 
  get_density(const typename Mesher_Gt::Point_3 &p) const {
    // NGHK: Make adaptive
    return 1;
  }

private:
  // Used to optimize the point location in TMC:
  mutable TMC_Cell_handle last_ch;

  void construct_bounding_box(Regular &regular);

  Gt &gt;
  RT shrink;
  Triangulated_mixed_complex _tmc;
  bool verbose;
  Sphere_3 _bounding_sphere;
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

    // Set the bounding sphere for the Delaunay mesher
    _bounding_sphere = Sphere_3(mid, dr*dr+1);
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

#endif // CGAL_SKIN_SURFACE_3_H
