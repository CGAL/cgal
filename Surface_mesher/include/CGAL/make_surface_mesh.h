// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_MAKE_SURFACE_MESH_H
#define CGAL_MAKE_SURFACE_MESH_H

#include <CGAL/Surface_mesher/Surface_mesher.h>
#include <CGAL/Surface_mesher/Surface_mesher_manifold.h>
#include <CGAL/Surface_mesher/Surface_mesher_edges_level.h>
#include <CGAL/Surface_mesher/Surface_mesher_visitor.h>
#include <CGAL/Surface_mesher/Surface_mesher_edges_level_visitor.h>
#include <CGAL/Surface_mesh_traits_generator_3.h>
#include <CGAL/Surface_mesh_complex_2_in_triangulation_3.h>

#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Surface_mesh_vertex_base_3.h>

#include <CGAL/Surface_mesh_default_criteria_3.h>

#include <CGAL/Surface_mesher/Has_edges.h>

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
                       const Surface& surface,
                       const Criteria& criteria,
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
                       const typename SurfaceMeshTraits_3::Surface_3& surface,
		       const SurfaceMeshTraits_3& surface_mesh_traits,
                       const Criteria& criteria,
                       Tag,
                       int initial_number_of_points = 20)
{
  typedef typename Make_surface_mesh_helper<
    C2T3,
    SurfaceMeshTraits_3,
    Criteria,
    Tag>::Mesher_base Mesher_base;

#ifdef CGAL_SURFACE_MESHER_VERBOSE
#error
  typedef Surface_mesher::Surface_mesher<
    Mesher_base,
    typename Surface_mesher::details::Facet_generator<Mesher_base>::type,
    Null_mesher_level,
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

/** For surfaces with curve edges. Implicit SurfaceMeshTraits_3. */
template <typename C2T3,
	  typename Surface,
	  typename FacetsCriteria,
	  typename EdgesCriteria,
	  typename Tag>
void make_piecewise_smooth_surface_mesh(C2T3& c2t3,
					const Surface& surface,
					const FacetsCriteria& facets_criteria,
					const EdgesCriteria& edges_criteria,
					Tag tag,
					int initial_number_of_points = 20)
{
  typedef typename Surface_mesh_traits_generator_3<Surface>::type Traits;

  make_piecewise_smooth_surface_mesh(c2t3, surface, Traits(), 
				     facets_criteria, edges_criteria,
				     tag,
				     initial_number_of_points);  
}

/** For surfaces with curve edges. Explicit SurfaceMeshTraits_3. */
template <typename C2T3,
	  typename SurfaceMeshTraits_3,
	  typename FacetsCriteria,
	  typename EdgesCriteria,
          typename Tag>
void 
make_piecewise_smooth_surface_mesh(C2T3& c2t3,
				   const typename SurfaceMeshTraits_3::Surface_3& surface,
				   const SurfaceMeshTraits_3& surface_mesh_traits,
				   const FacetsCriteria& facets_criteria,
				   const EdgesCriteria& edges_criteria,
				   Tag,
				   int initial_number_of_points = 20)
{
  typedef typename Make_surface_mesh_helper<
    C2T3,
    SurfaceMeshTraits_3,
    FacetsCriteria,
    Tag>::Mesher_base Facets_mesher_base;

  typedef Surface_mesher::Surface_mesher_edges_level_base<
    C2T3,
    typename SurfaceMeshTraits_3::Surface_3,
    SurfaceMeshTraits_3,
    EdgesCriteria> Edges_mesher_base;

#ifdef CGAL_SURFACE_MESHER_VERBOSE
  typedef Surface_mesher::Surface_mesher<
    Edges_mesher_base,
    typename Surface_mesher::details::Edge_generator<Edges_mesher_base>::type,
    Null_mesher_level,
    Surface_mesher::VERBOSE> Edges_level;
#else
  typedef Surface_mesher::Surface_mesher<
    Edges_mesher_base,
    typename Surface_mesher::details::Edge_generator<Edges_mesher_base>::type
  > Edges_level;
#endif

#ifdef CGAL_SURFACE_MESHER_VERBOSE
  typedef Surface_mesher::Surface_mesher<
    Facets_mesher_base,
    typename Surface_mesher::details::Facet_generator<Facets_mesher_base>::type,
    Edges_level,
    Surface_mesher::VERBOSE> Facets_level;
#else
  typedef Surface_mesher::Surface_mesher<
    Facets_mesher_base,
    typename Surface_mesher::details::Facet_generator<Facets_mesher_base>::type,
    Edges_level> Facets_level;
#endif

  typename SurfaceMeshTraits_3::Construct_initial_points get_initial_points =
    surface_mesh_traits.construct_initial_points_object();

  typedef Surface_mesher::Visitor<typename C2T3::Triangulation, 
    Facets_level, Null_mesh_visitor> Edges_level_visitor;
  typedef Surface_mesher::Edges_level_visitor<typename C2T3::Triangulation, 
    Edges_level, Edges_level_visitor> Facets_and_edges_visitor;

  get_initial_points(surface,
                     CGAL::inserter(c2t3.triangulation()),
                     initial_number_of_points);
  Edges_level edges_level(c2t3, surface, surface_mesh_traits, edges_criteria);
  Facets_level facets_level(c2t3, surface, surface_mesh_traits, facets_criteria,
			    edges_level);

  Null_mesh_visitor null_mesh_visitor;
  Edges_level_visitor edges_level_visitor(&facets_level, &null_mesh_visitor);
  Facets_and_edges_visitor facets_and_edges_visitor(&edges_level, &edges_level_visitor);

  edges_level.init();
  edges_level.refine_mesh();
  facets_level.scan_triangulation();
  facets_level.refine_mesh(facets_and_edges_visitor);
}

} // end namespace CGAL

#endif // CGAL_MAKE_SURFACE_MESH_H
