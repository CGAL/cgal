#ifndef CGAL_DEMO_MESH_3_C3T3_TYPE_H
#define CGAL_DEMO_MESH_3_C3T3_TYPE_H

#include "Polyhedron_type.h"
#include "Image_type.h"
#include "implicit_functions/Implicit_function_interface.h"

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
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

namespace CGAL {

// A specialisation of Triangle_accessor_3 which fits our Polyhedron type
template <typename K>
class Triangle_accessor_3<Polyhedron, K>
{
public:
  typedef typename K::Triangle_3                    Triangle_3;
  typedef Polyhedron::Facet_const_iterator Triangle_iterator;
  typedef Polyhedron::Facet_const_handle   Triangle_handle;
  
  Triangle_accessor_3() { }
  
  Triangle_iterator triangles_begin(const Polyhedron& p) const
  {
    return p.facets_begin();
  }
  
  Triangle_iterator triangles_end(const Polyhedron& p) const
  {
    return p.facets_end();
  }
  
  Triangle_3 triangle(const Triangle_handle& handle) const
  {
    typedef typename K::Point_3 Point;
    const Point& a = handle->halfedge()->vertex()->point();
    const Point& b = handle->halfedge()->next()->vertex()->point();
    const Point& c = handle->halfedge()->next()->next()->vertex()->point();
    return Triangle_3(a,b,c);
  }
};

}

typedef CGAL::Mesh_3::Robust_intersection_traits_3<Kernel>              RKernel;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron,RKernel>              Polyhedral_mesh_domain;
typedef CGAL::Labeled_image_mesh_domain_3<Image,Kernel>                 Image_mesh_domain;
typedef Wrapper<Kernel>                                                 Function_wrapper;
typedef CGAL::Mesh_3::Labeled_mesh_domain_3<Function_wrapper, Kernel>   Function_mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Polyhedral_mesh_domain>::type Tr;

// 3D complex
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

#endif // CGAL_DEMO_MESH_3_C3T3_TYPE_H
