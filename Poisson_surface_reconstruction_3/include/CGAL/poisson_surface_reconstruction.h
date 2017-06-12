// Copyright (c) 2017  GeometryFactory (France)
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_POISSON_SURFACE_RECONSTRUCTION_H
#define CGAL_POISSON_SURFACE_RECONSTRUCTION_H

#include <CGAL/license/Poisson_surface_reconstruction_3.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/property_map.h>

namespace CGAL {

  
  /*!
    \ingroup PkgPoissonSurfaceReconstruction

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
    point set (using `compute_average_spacing()` for example). Higher
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
#if defined(DOXYGEN_RUNNING) || !defined(CGAL_CFG_NO_CPP0X_DEFAULT_TEMPLATE_ARGUMENTS_FOR_FUNCTION_TEMPLATES)
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
#else
  template <typename PointInputIterator,
            typename PointMap,
            typename NormalMap,
            typename PolygonMesh>
  bool
  poisson_surface_reconstruction_delaunay (PointInputIterator begin,
                                           PointInputIterator end,
                                           PointMap point_map,
                                           NormalMap normal_map,
                                           PolygonMesh& output_mesh,
                                           double spacing,
                                           double sm_angle = 20.0,
                                           double sm_radius = 30.0,
                                           double sm_distance = 0.375)
  {
    return poisson_surface_reconstruction_delaunay (begin, end, point_map, normal_map, output_mesh,
                                                    spacing, sm_angle, sm_radius, sm_distance,
                                                    CGAL::Manifold_with_boundary_tag());
  }

  template <typename PointInputIterator,
            typename PointMap,
            typename NormalMap,
            typename PolygonMesh,
            typename Tag>
  bool
  poisson_surface_reconstruction_delaunay (PointInputIterator begin,
                                           PointInputIterator end,
                                           PointMap point_map,
                                           NormalMap normal_map,
                                           PolygonMesh& output_mesh,
                                           double spacing,
                                           double sm_angle = 20.0,
                                           double sm_radius = 30.0,
                                           double sm_distance = 0.375)
  {
    return poisson_surface_reconstruction_delaunay (begin, end, point_map, normal_map, output_mesh,
                                                    spacing, sm_angle, sm_radius, sm_distance,
                                                    Tag());
  }

  template <typename PointInputIterator,
            typename PointMap,
            typename NormalMap,
            typename PolygonMesh,
            typename Tag>
  bool
  poisson_surface_reconstruction_delaunay (PointInputIterator begin,
                                           PointInputIterator end,
                                           PointMap point_map,
                                           NormalMap normal_map,
                                           PolygonMesh& output_mesh,
                                           double spacing,
                                           double sm_angle,
                                           double sm_radius,
                                           double sm_distance,
                                           Tag tag)
#endif
  {
    typedef typename boost::property_traits<PointMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    typedef typename Kernel::Sphere_3 Sphere;
    
    typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
    typedef CGAL::Surface_mesh_default_triangulation_3 STr;
    typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
    typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;
    
    Poisson_reconstruction_function function(begin, end, point_map, normal_map);
    if ( ! function.compute_implicit_function() ) 
      return false;

    Point inner_point = function.get_inner_point();
    Sphere bsphere = function.bounding_sphere();
    double radius = std::sqrt(bsphere.squared_radius());

    double sm_sphere_radius = 5.0 * radius;
    double sm_dichotomy_error = sm_distance * spacing / 1000.0;
    
    Surface_3 surface(function,
                      Sphere (inner_point, sm_sphere_radius * sm_sphere_radius),
                      sm_dichotomy_error / sm_sphere_radius);

    CGAL::Surface_mesh_default_criteria_3<STr> criteria (sm_angle,
                                                         sm_radius * spacing,
                                                         sm_distance * spacing);

    STr tr;
    C2t3 c2t3(tr);
    
    CGAL::make_surface_mesh(c2t3,
                            surface,
                            criteria,
                            tag);

    if(tr.number_of_vertices() == 0)
      return false;

    CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);

    return true;
  }


}


#endif // CGAL_POISSON_SURFACE_RECONSTRUCTION_H
