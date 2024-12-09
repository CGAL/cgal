// Copyright (c) 2017, 2024  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot, Jane Tournois

#ifndef CGAL_POISSON_SURFACE_RECONSTRUCTION_H
#define CGAL_POISSON_SURFACE_RECONSTRUCTION_H

#include <CGAL/license/Poisson_surface_reconstruction_3.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/property_map.h>

namespace CGAL {


  /*!
    \ingroup PkgPoissonSurfaceReconstruction3Ref

    Performs surface reconstruction as follows:

    - compute the Poisson implicit function, through a conjugate
      gradient solver, represented as a piecewise linear function
      stored on a 3D Delaunay mesh generated via Delaunay refinement
    - meshes the function with a user-defined precision using another
      round of Delaunay refinement: it contours the isosurface
      corresponding to the isovalue of the median of the function
      values at the input points
    - outputs the result in a polygon mesh

    This function relies mainly on the size parameter `spacing`. A
    reasonable solution is to use the average spacing of the input
    point set (using `compute_average_spacing()` for example). Smaller
    values increase the precision of the output mesh at the cost of
    higher computation time.

    Parameters `sm_angle`, `sm_radius` and `sm_distance` work
    similarly to the parameters of `SurfaceMeshFacetsCriteria_3`. The
    latest two are defined with respect to `spacing`.

    \tparam PointInputIterator is a model of `InputIterator`.

    \tparam PointMap is a model of `ReadablePropertyMap` with value
    type `Point_3<Kernel>`.

    \tparam NormalMap is a model of `ReadablePropertyMap` with value
    type `Vector_3<Kernel>`.

    \tparam PolygonMesh a model of `MutableFaceGraph` with an internal
    point property map.

    \tparam Tag is a tag whose type affects the behavior of the
    meshing algorithm (see `make_surface_mesh()`).

    \param begin iterator on the first point of the sequence.
    \param end past the end iterator of the point sequence.
    \param point_map property map: value_type of `InputIterator` -> Point_3.
    \param normal_map property map: value_type of `InputIterator` -> Vector_3.
    \param output_mesh where the reconstruction is stored.
    \param spacing size parameter.
    \param sm_angle bound for the minimum facet angle in degrees.
    \param sm_radius bound for the radius of the surface Delaunay balls (relatively to the `average_spacing`).
    \param sm_distance bound for the center-center distances (relatively to the `average_spacing`).
    \param tag surface mesher tag.
    \return `true` if reconstruction succeeded, `false` otherwise.
  */
  template <typename PointInputIterator,
            typename PointMap,
            typename NormalMap,
            typename PolygonMesh,
            typename Tag = CGAL::Manifold_with_boundary_tag>
  bool
  poisson_surface_reconstruction_delaunay (PointInputIterator begin,
                                           PointInputIterator end,
                                           PointMap point_map,
                                           NormalMap normal_map,
                                           PolygonMesh& output_mesh,
                                           double spacing,
                                           double sm_angle = 20.0,
                                           double sm_radius = 30.0,
                                           double sm_distance = 0.375,
                                           Tag tag = Tag())
  {
    typedef typename boost::property_traits<PointMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    typedef typename Kernel::Sphere_3 Sphere;
    typedef typename Kernel::FT FT;

    typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
    typedef CGAL::Labeled_mesh_domain_3<Kernel> Mesh_domain;
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Sequential_tag>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

    Poisson_reconstruction_function function(begin, end, point_map, normal_map);
    if ( ! function.compute_implicit_function() )
      return false;

    Point inner_point = function.get_inner_point();
    Sphere bsphere = function.bounding_sphere();
    FT radius = CGAL::approximate_sqrt(bsphere.squared_radius());

    FT sm_sphere_radius = 5.0 * radius;
    FT sm_dichotomy_error = sm_distance * spacing / 1000.0;

    Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain(function, Sphere(inner_point, sm_sphere_radius),
      CGAL::parameters::relative_error_bound(sm_dichotomy_error / sm_sphere_radius));

    Mesh_criteria criteria(CGAL::parameters::facet_angle = sm_angle,
                           CGAL::parameters::facet_size = sm_radius*spacing,
                           CGAL::parameters::facet_distance = sm_distance*spacing);


    auto turn_tag_into_mesh_3_manifold_option = [](Tag) {
      if constexpr (std::is_same_v<Tag, CGAL::Manifold_with_boundary_tag>)
        return CGAL::parameters::manifold_with_boundary();
      else if constexpr (std::is_same_v<Tag, CGAL::Manifold_tag>)
        return CGAL::parameters::manifold();
      else
        return CGAL::parameters::non_manifold();
    };

    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        turn_tag_into_mesh_3_manifold_option(tag)
                                        .no_exude().no_perturb()
                                        .manifold_with_boundary());

    const auto& tr = c3t3.triangulation();

    if(tr.number_of_vertices() == 0)
      return false;

    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, output_mesh);

    return true;
  }


}


#endif // CGAL_POISSON_SURFACE_RECONSTRUCTION_H
