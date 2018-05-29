#ifndef CGAL_DEMO_MESH_3_C3T3_TYPE_H
#define CGAL_DEMO_MESH_3_C3T3_TYPE_H

// include this to get #define BOOST_PARAMETER_MAX_ARITY 12
// as otherwise it gets set via inclusion of Polyhedron_3.h
#include <CGAL/boost/parameter.h>
#include <CGAL/Default.h>

#include "Polyhedron_type.h"
#include "SMesh_type.h"
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
#include "Image_type.h"
#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include "Polyhedron_demo_mesh_3_labeled_mesh_domain_3.h"
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
#include "implicit_functions/Implicit_function_interface.h"
#include "Polyhedron_demo_mesh_3_labeled_mesh_domain_3.h"
#endif

#include <CGAL/Mesh_triangulation_3.h>
#include "Scene_surface_mesh_item.h"
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>

#include <CGAL/tags.h>

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
template <typename K>
struct Wrapper
{
  typedef int return_type;
  typedef typename K::Point_3 Point_3;
  
  Wrapper(const Implicit_function_interface& f) : f_(f) {}
  return_type operator()(const Point_3& p, const bool=true) const
  {
    return (f_(p.x(),p.y(),p.z()) < 0) ? 1 : 0;
  }
  
private:
  const Implicit_function_interface& f_;
};
#endif

typedef CGAL::Polyhedral_mesh_domain_with_features_3<
          Kernel, Polyhedron, CGAL::Default, CGAL::Tag_true> Polyhedral_mesh_domain;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<
          Kernel, SMesh, CGAL::Default, int> Polyhedral_mesh_domain_sm;
// The last `Tag_true` says the Patch_id type will be int, and not pair<int, int>

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
typedef CGAL::Labeled_image_mesh_domain_3<Image,Kernel>           Image_domain1;
typedef CGAL::Polyhedron_demo_labeled_mesh_domain_3<Image_domain1> Image_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Image_domain> Image_mesh_domain;
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
typedef Wrapper<Kernel>                                              Function_wrapper;
typedef CGAL::Labeled_mesh_domain_3<Function_wrapper, Kernel>        Function_domain;
typedef CGAL::Polyhedron_demo_labeled_mesh_domain_3<Function_domain> Function_mesh_domain;
#endif

//Robust_cc_geom_traits
typedef CGAL::Kernel_traits<Polyhedral_mesh_domain>::Kernel
    Robust_intersections_traits;
typedef CGAL::details::Mesh_geom_traits_generator<Robust_intersections_traits>::type
    Robust_K;

// Triangulation
typedef CGAL::Compact_mesh_cell_base_3<Robust_K, Polyhedral_mesh_domain>    Cell_base;
typedef CGAL::Triangulation_cell_base_with_info_3<int, Robust_K, Cell_base> Cell_base_with_info;

#ifdef CGAL_CONCURRENT_MESH_3
  typedef CGAL::Mesh_triangulation_3<Polyhedral_mesh_domain, 
                                     Robust_intersections_traits,
                                     CGAL::Parallel_tag,
                                     CGAL::Default,
                                     Cell_base_with_info>::type Tr;
#else
  typedef CGAL::Mesh_triangulation_3<Polyhedral_mesh_domain,
                                     Robust_intersections_traits,
                                     CGAL::Sequential_tag,
                                     CGAL::Default,
                                     Cell_base_with_info>::type Tr;
#endif

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

typedef Tr::Geom_traits Geom_traits;

// 3D complex

#endif // CGAL_DEMO_MESH_3_C3T3_TYPE_H
