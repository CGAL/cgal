// Copyright (c) 2005 RuG (Netherlands)
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
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_SKIN_SURFACE_TRAITS_3_H
#define CGAL_SKIN_SURFACE_TRAITS_3_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// Contains the weighted converter:
#include <CGAL/Regular_triangulation_filtered_traits_3.h>
//#include <CGAL/Regular_triangulation_3.h>

//#include <CGAL/Skin_surface_simplicial_complex_3.h>

//#include <CGAL/Polyhedron_3.h>
//#include <CGAL/Skin_surface_polyhedral_items_3.h>

CGAL_BEGIN_NAMESPACE 

template <class Kernel>
class Construct_anchor_point_3 {
public:
  typedef typename Kernel::RT            RT;
  typedef typename Kernel::Point_3       Point;
  
  Construct_anchor_point_3(RT const &shrink) : shrink(shrink) {}
  
  Point operator()(const Point &center_del, const Point &center_vor) {
    return center_del + shrink*(center_vor - center_del);
  }
  RT const shrink;
};

// The following four converters convert geometric objects between the
//   datastructures mentioned above:
// - RegToSimplConverter
// - RegToMeshConverter 
// - SimplToMeshConverter 
// - MeshToSimplConverter 

template <
  class RegularKernel = Exact_predicates_inexact_constructions_kernel,
  class TriangulatedMixedComplexKernel =
  //Exact_predicates_inexact_constructions_kernel,
  Exact_predicates_exact_constructions_kernel,
  class PolyhedronKernel = Simple_cartesian<double> >  
class Skin_surface_traits_3 {
public:
  typedef Skin_surface_traits_3<RegularKernel, TriangulatedMixedComplexKernel,
				PolyhedronKernel>                         Self;
  // Kernel definitions and converters
  typedef RegularKernel                      Regular_kernel;
  typedef Regular_triangulation_euclidean_traits_3<Regular_kernel>
                                             Regular_traits;
  typedef TriangulatedMixedComplexKernel     Triangulated_mixed_complex_traits;
  typedef PolyhedronKernel                   Polyhedron_traits;

  typedef typename Regular_kernel::RT        Regular_RT;

  typedef Regular_triangulation_euclidean_traits_3<
           Triangulated_mixed_complex_traits> Weighted_triangulated_mixed_complex_traits;

  typedef typename Weighted_triangulated_mixed_complex_traits::
    Construct_weighted_circumcenter_3         Construct_weighted_circumcenter_3;

  typedef Weighted_converter_3 <
    Cartesian_converter<Regular_kernel,Triangulated_mixed_complex_traits> >
                                              R2T_converter;
  typedef Weighted_converter_3 <
    Cartesian_converter<Regular_kernel, Polyhedron_traits> >
                                              R2P_converter;
  typedef Weighted_converter_3 <
    Cartesian_converter<Triangulated_mixed_complex_traits,Polyhedron_traits> >
                                              T2P_converter; 
  typedef Weighted_converter_3 <
    Cartesian_converter<Polyhedron_traits,Triangulated_mixed_complex_traits> >
                                              P2T_converter; 

  Skin_surface_traits_3(Regular_RT shrink = .5) :
    shrink(shrink) {}
  
  // Function objects:
  R2T_converter r2t_converter_object() const {
    return R2T_converter();
  }
  R2P_converter r2p_converter_object() const {
    return R2P_converter();
  }
  T2P_converter t2p_converter_object() const {
    return T2P_converter();
  }
  P2T_converter p2t_converter_object() const {
    return P2T_converter();
  }

//   template <class Kernel>
//   Construct_weighted_circumcenter_3<
//     Regular_triangulation_euclidean_traits_3<Kernel> >
//   construct_weighted_circumcenter_3_object() const {
//     return Regular_triangulation_euclidean_traits_3<Kernel>().
//       construct_weighted_circumcenter_3_object();
//   }

  template <class Kernel>
  In_smallest_orthogonal_sphere_3<
    Regular_triangulation_euclidean_traits_3<Kernel> >
  in_smallest_orthogonal_sphere_3_object() const {
    return Regular_triangulation_euclidean_traits_3<Kernel>().
      in_smallest_orthogonal_sphere_3_object();
  }

  Construct_anchor_point_3<Triangulated_mixed_complex_traits>
  construct_anchor_point_3_object() const {
    return Construct_anchor_point_3<Triangulated_mixed_complex_traits>(
      r2t_converter_object()(shrink) );
  }
  Regular_RT shrink_factor() const {
    return shrink;
  }

  Regular_RT shrink;
};

CGAL_END_NAMESPACE

#endif // CGAL_SKIN_SURFACE_TRAITS_3_H
