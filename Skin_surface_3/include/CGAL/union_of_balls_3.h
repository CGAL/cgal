// NGHK: NOT IN USE FOR NOW




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

#include <CGAL/Bbox_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulated_mixed_complex_3.h>
#include <CGAL/triangulate_voronoi_diagram_3.h>
#include <CGAL/Marching_tetrahedra_traits_skin_surface_3.h>
#include <CGAL/Marching_tetrahedra_observer_default_3.h>
#include <CGAL/marching_tetrahedra_3.h>
#include <CGAL/skin_surface_3.h>

CGAL_BEGIN_NAMESPACE

template <class InputIterator, class Polyhedron_3, class SkinSurfaceTraits_3>
void union_of_balls_3(InputIterator first, InputIterator last,
  Polyhedron_3 &polyhedron, const SkinSurfaceTraits_3 &skin_surface_traits,
  bool verbose = false) {
  std::cout << "NGHK: NOT YET WORKING";
  if (first == last) {
    std::cout << " No input balls" << std::endl;
    return;
  }

  // Types
  typedef SkinSurfaceTraits_3                              Skin_surface_traits;
  typedef typename Skin_surface_traits::Regular_traits     Regular_traits;
  typedef typename Regular_traits::Bare_point              Reg_point;
  typedef typename Regular_traits::Weighted_point          Reg_weighted_point;

  typedef Regular_triangulation_3<Regular_traits> Regular;
  typedef Triangulated_mixed_complex_3<SkinSurfaceTraits_3>
                                                  Triangulated_mixed_complex;
  typedef Marching_tetrahedra_traits_skin_surface_3<
    Triangulated_mixed_complex,
    Polyhedron_3,
    typename SkinSurfaceTraits_3::T2P_converter>  Marching_tetrahedra_traits;
  typedef Marching_tetrahedra_observer_default_3<
    Triangulated_mixed_complex, Polyhedron_3>     Marching_tetrahedra_observer;
    
  // Code
  Regular regular;
  Triangulated_mixed_complex triangulated_mixed_complex;

  while (first != last) {
    regular.insert((*first));
    first++;
  }

  skin_surface_construct_bounding_box_3(regular,skin_surface_traits);
  
  if (verbose) {
    std::cerr << "Triangulation ready" << std::endl;
  }

  // Construct the triangulated mixed complex:
  triangulate_voronoi_diagram_3(
    regular, triangulated_mixed_complex, skin_surface_traits);

  CGAL_assertion(triangulated_mixed_complex.is_valid());
  if (verbose) {
    std::cerr << "Triangulated mixed complex ready" << std::endl;
  }

  // Extract the coarse mesh using marching_tetrahedra
  Marching_tetrahedra_traits marching_traits;
  marching_tetrahedra_3(
    triangulated_mixed_complex, polyhedron, marching_traits);

  if (verbose) {
    std::cerr << "Mesh ready" << std::endl;
  }
  
}


CGAL_END_NAMESPACE

#endif // CGAL_UNION_OF_BALLS_3_H
