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

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/property_map.h>
#include <CGAL/compute_average_spacing.h>

namespace CGAL {

  
  template <typename ConcurrencyTag,
            typename PointInputIterator,
            typename PointMap,
            typename NormalMap,
            typename PolygonMesh>
  bool
  poisson_surface_reconstruction(PointInputIterator begin,
                                 PointInputIterator end,
                                 PointMap point_map,
                                 NormalMap normal_map,
                                 PolygonMesh& output_mesh,
                                 double sm_angle = 20.0,
                                 double sm_radius = 30.0,
                                 double sm_distance = 0.375)
  {
    typedef typename boost::property_traits<PointMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Sphere_3 Sphere;
    
    typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
    typedef CGAL::Surface_mesh_default_triangulation_3 STr;
    typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
    typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;
    
    Poisson_reconstruction_function function(begin, end, point_map, normal_map);
    if ( ! function.compute_implicit_function() ) 
      return false;

    FT average_spacing = CGAL::compute_average_spacing<ConcurrencyTag>
      (begin, end, point_map, 6);

    Point inner_point = function.get_inner_point();
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());

    FT sm_sphere_radius = 5.0 * radius;
    FT sm_dichotomy_error = sm_distance * average_spacing / 1000.0;
    
    Surface_3 surface(function,
                      Sphere (inner_point, sm_sphere_radius * sm_sphere_radius),
                      sm_dichotomy_error / sm_sphere_radius);

    CGAL::Surface_mesh_default_criteria_3<STr> criteria (sm_angle,
                                                         sm_radius * average_spacing,
                                                         sm_distance * average_spacing);

    STr tr;
    C2t3 c2t3(tr);
    
    CGAL::make_surface_mesh(c2t3,
                            surface,
                            criteria,
                            CGAL::Manifold_with_boundary_tag());

    if(tr.number_of_vertices() == 0)
      return false;

    CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);

    return true;
  }


}


#endif // CGAL_POISSON_SURFACE_RECONSTRUCTION_H
