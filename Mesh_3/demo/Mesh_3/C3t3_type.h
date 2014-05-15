#ifndef CGAL_DEMO_MESH_3_C3T3_TYPE_H
#define CGAL_DEMO_MESH_3_C3T3_TYPE_H

#include "Polyhedron_type.h"
#include "Image_type.h"
#include "implicit_functions/Implicit_function_interface.h"

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Labeled_image_mesh_domain_3.h>

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


typedef CGAL::Triangle_accessor_3<Polyhedron, Kernel> T_accessor;

typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel,
                                                     Polyhedron,
                                                     T_accessor,
                                                     CGAL::Tag_false>
                                                                        Polyhedral_mesh_domain;
// The last `Tag_false` says the Patch_id type will be pair<int,int>, like
// the `Image_mesh_domain`, and `Function_mesh_domain`.

typedef CGAL::Labeled_image_mesh_domain_3<Image,Kernel>                 Image_mesh_domain;
typedef Wrapper<Kernel>                                                 Function_wrapper;
typedef CGAL::Mesh_3::Labeled_mesh_domain_3<Function_wrapper, Kernel>   Function_mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Polyhedral_mesh_domain>::type Tr;

// 3D complex
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

#endif // CGAL_DEMO_MESH_3_C3T3_TYPE_H
