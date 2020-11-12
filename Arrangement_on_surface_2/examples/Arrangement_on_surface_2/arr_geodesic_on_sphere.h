#ifndef ARR_GEODESIC_ON_SPHERE_H
#define ARR_GEODESIC_ON_SPHERE_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel  Kernel;
typedef Kernel::FT                                         Number_type;

typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>  Geom_traits;
typedef Geom_traits::Point_2                               Point;
typedef Geom_traits::X_monotone_curve_2                    X_monotone_curve;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits> Topol_traits;
typedef CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>
                                                           Arrangement;

typedef Arrangement::Vertex_handle                         Vertex_handle;
typedef Arrangement::Halfedge_handle                       Halfedge_handle;
typedef Arrangement::Face_handle                           Face_handle;
typedef Arrangement::Vertex_const_handle                   Vertex_const_handle;
typedef Arrangement::Halfedge_const_handle                 Halfedge_const_handle;
typedef Arrangement::Face_const_handle                     Face_const_handle;

#endif
