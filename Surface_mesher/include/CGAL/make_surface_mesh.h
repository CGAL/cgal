// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
// $Revision$ $Date$
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MAKE_SURFACE_MESH_H
#define CGAL_MAKE_SURFACE_MESH_H

#include <CGAL/Surface_mesher/Surface_mesher.h>
#include <CGAL/Surface_mesher/Surface_mesher_manifold.h>
#include <CGAL/Surface_mesh_traits_generator_3.h>
#include <CGAL/Surface_mesh_complex_2_in_triangulation_3.h>

#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Surface_mesh_vertex_base_3.h>

#include <CGAL/Surface_mesh_default_criteria_3.h>

#include <CGAL/iterator.h> // CGAL::inserter()

namespace CGAL {

  struct Non_manifold_tag {};
  struct Manifold_tag {};
  struct Manifold_with_boundary_tag {};

  template <
    typename C2T3,
    typename SurfaceMeshTraits_3,
    typename Criteria,
    typename Tag // generic version: generates a compile time error
  >
  struct Make_surface_mesh_helper {
    template <typename T>
    struct Tag_does_not_exist_error {};
    typedef Tag_does_not_exist_error<Tag> Mesher_base;
  };

  template <
    typename C2T3,
    typename SurfaceMeshTraits_3,
    typename Criteria
  >
  struct Make_surface_mesh_helper<
    C2T3,
    SurfaceMeshTraits_3,
    Criteria,
    Non_manifold_tag> // Non_manifold_tag partial specialization
  {
    typedef Surface_mesher::Surface_mesher_base<
      C2T3,
      typename SurfaceMeshTraits_3::Surface_3,
      SurfaceMeshTraits_3,
      Criteria> Mesher_base;
  };

  template <
    typename C2T3,
    typename SurfaceMeshTraits_3,
    typename Criteria
  >
  struct Make_surface_mesh_helper<
    C2T3,
    SurfaceMeshTraits_3,
    Criteria,
    Manifold_with_boundary_tag> // Manifold_with_boundary_tag partial
                                // specialization
  {
    typedef Surface_mesher::Surface_mesher_regular_edges_base<
      C2T3,
      typename SurfaceMeshTraits_3::Surface_3,
      SurfaceMeshTraits_3,
      Criteria,
      true> /* true means "with boundary"*/ Regular_edge_base;

    typedef Surface_mesher::Surface_mesher_manifold_base<
      C2T3,
      typename SurfaceMeshTraits_3::Surface_3,
      SurfaceMeshTraits_3,
      Criteria,
      Regular_edge_base
      > Mesher_base;
  };

  template <
    typename C2T3,
    typename SurfaceMeshTraits_3,
    typename Criteria
  >
  struct Make_surface_mesh_helper<
    C2T3,
    SurfaceMeshTraits_3,
    Criteria,
    Manifold_tag> // Manifold_tag partial specialization
  {
    typedef Surface_mesher::Surface_mesher_regular_edges_base<
      C2T3,
      typename SurfaceMeshTraits_3::Surface_3,
      SurfaceMeshTraits_3,
      Criteria> Regular_edge_without_boundary_base;

    typedef Surface_mesher::Surface_mesher_manifold_base<
      C2T3,
      typename SurfaceMeshTraits_3::Surface_3,
      SurfaceMeshTraits_3,
      Criteria,
      Regular_edge_without_boundary_base
      > Mesher_base;
  };

template <typename C2T3,
	  typename Surface,
	  typename Criteria,
	  typename Tag>
void make_surface_mesh(C2T3& c2t3,
                       Surface& surface,
                       Criteria criteria,
                       Tag tag,
                       int initial_number_of_points = 20)  // TODO: document
                                                           // this parameter
{
  typedef typename Surface_mesh_traits_generator_3<Surface>::type Traits;

  make_surface_mesh(c2t3, surface, Traits(), criteria, tag,
                    initial_number_of_points);  
}

template <typename C2T3,
	  typename SurfaceMeshTraits_3,
          typename Criteria,
          typename Tag>
void make_surface_mesh(C2T3& c2t3,
                       typename SurfaceMeshTraits_3::Surface_3& surface,
		       SurfaceMeshTraits_3 surface_mesh_traits,
                       Criteria criteria,
                       Tag,
                       int initial_number_of_points = 20)
{
  typedef typename Make_surface_mesh_helper<
    C2T3,
    SurfaceMeshTraits_3,
    Criteria,
    Tag>::Mesher_base Mesher_base;

#ifdef CGAL_SURFACE_MESHER_VERBOSE
  typedef Surface_mesher::Surface_mesher<
    Mesher_base,
    Surface_mesher::VERBOSE> Mesher;
#else
  typedef Surface_mesher::Surface_mesher<Mesher_base> Mesher;
#endif

  typename SurfaceMeshTraits_3::Construct_initial_points get_initial_points =
    surface_mesh_traits.construct_initial_points_object();

  get_initial_points(surface,
                     CGAL::inserter(c2t3.triangulation()),
                     initial_number_of_points);
  Mesher mesher(c2t3, surface, surface_mesh_traits, criteria);
  mesher.refine_mesh();
  // TODO initial, then refine()
}

} // end namespace CGAL

#endif // CGAL_MAKE_SURFACE_MESH_H
