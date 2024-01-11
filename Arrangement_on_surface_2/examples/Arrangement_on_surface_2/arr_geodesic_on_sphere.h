#ifndef ARR_GEODESIC_ON_SPHERE_H
#define ARR_GEODESIC_ON_SPHERE_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Number_type = Kernel::FT;

using Geom_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
using Point = Geom_traits::Point_2;
using X_monotone_curve = Geom_traits::X_monotone_curve_2;
using Topol_traits = CGAL::Arr_spherical_topology_traits_2<Geom_traits>;
using Arrangement = CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>;

using Vertex_handle = Arrangement::Vertex_handle;
using Halfedge_handle = Arrangement::Halfedge_handle;
using Face_handle = Arrangement::Face_handle;
using Vertex_const_handle = Arrangement::Vertex_const_handle;
using Halfedge_const_handle = Arrangement::Halfedge_const_handle;
using Face_const_handle = Arrangement::Face_const_handle;

#endif
