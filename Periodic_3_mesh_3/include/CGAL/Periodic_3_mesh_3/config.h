#ifndef CGAL_PERIODIC_3_MESH_3_CONFIG_H
#define CGAL_PERIODIC_3_MESH_3_CONFIG_H

// Avoid optimisations of Mesh_3
#define CGAL_NO_STRUCTURAL_FILTERING
#ifdef CGAL_MESH_3_SIZING_FIELD_INEXACT_LOCATE
  #undef CGAL_MESH_3_SIZING_FIELD_INEXACT_LOCATE
#endif

// to be sure that we include it before Mesh_3's version
#include <CGAL/Periodic_3_mesh_3/Protect_edges_sizing_field.h>

//#warning "Structural filtering is disabled because this feature isn't compatible with some periodic traits!"
/*
 * It's a temporary solution.
 * We avoid the structural filtering optimization, because it isn't
 * available with some periodic triangulation traits. (cf. Periodic_3_triangulation_remove_traits_3.h)
 * Indeed, this feature needs the Triangulation_3::inexact_orientation() method.
 *
 * A better solution would be to make Triangulation_3 compatible 'any' traits.
 * By this way, Triangulation_3 would require the Bare_point concept which would use in
 * inexact_orientation(). So, we'd work with a Bare_point instead of a Point in this method.
 * In the 'default' Triangulation_3 traits, Bare_point would be equal to Point_3
 * (or const Point_3&  <- for avoiding useless copies).
*/

#endif // CGAL_PERIODIC_3_MESH_3_CONFIG_H
