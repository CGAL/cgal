#ifndef CGAL_TEST_EXTENDED_TOPOL_TRAITS_H
#define CGAL_TEST_EXTENDED_TOPOL_TRAITS_H

// ============================================================================
// Topology Traits includes:
// ============================================================================
#if TEST_TOPOL_TRAITS == PLANAR_BOUNDED_TOPOL_TRAITS
#include <CGAL/Arr_bounded_planar_topology_traits_2.h>

#elif TEST_TOPOL_TRAITS == PLANAR_UNBOUNDED_TOPOL_TRAITS
#include <CGAL/Arr_unb_planar_topology_traits_2.h>

#elif TEST_TOPOL_TRAITS == SPHERICAL_TOPOL_TRAITS
#include <CGAL/Arr_spherical_topology_traits_2.h>

#else
#error No topology traits (TOPOL TRAITS) specified!
#endif

// ============================================================================
// Topology Traits typedef:
// ============================================================================
#if TEST_TOPOL_TRAITS == PLANAR_BOUNDED_TOPOL_TRAITS
typedef CGAL::Arr_bounded_planar_topology_traits_2<Geom_traits, Dcel>
                                                                Topol_traits;
#define TOPOL_TRAITS_TYPE "Bounded Planar"

#elif TEST_TOPOL_TRAITS == PLANAR_UNBOUNDED_TOPOL_TRAITS
typedef CGAL::Arr_unb_planar_topology_traits_2<Geom_traits, Dcel>
                                                                Topol_traits;
#define TOPOL_TRAITS_TYPE "Unbounded Planar"

#elif TEST_TOPOL_TRAITS == SPHERICAL_TOPOL_TRAITS
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits, Dcel>
                                                                Topol_traits;
#define TOPOL_TRAITS_TYPE "Spherical"

#else
#error No topology traits (TOPOL TRAITS) specified!
#endif

#endif
