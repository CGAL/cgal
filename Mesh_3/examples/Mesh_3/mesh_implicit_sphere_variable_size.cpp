#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Implicit_mesh_domain_3<Function,K> Mesh_domain;

// Triangulation
#ifdef CGAL_CONCURRENT_MESH_3
  typedef CGAL::Mesh_triangulation_3<
    Mesh_domain,
    CGAL::Kernel_traits<Mesh_domain>::Kernel, // Same as sequential
    CGAL::Parallel_tag                        // Tag to activate parallelism
  >::type Tr;
#else
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// Sizing field
struct Spherical_sizing_field
{
  typedef ::FT FT;
  typedef Point Point_3;
  typedef Mesh_domain::Index Index;
  
  FT operator()(const Point_3& p, const int, const Index&) const
  {
    FT sq_d_to_origin = CGAL::squared_distance(p, Point(CGAL::ORIGIN));
    return CGAL::abs( CGAL::sqrt(sq_d_to_origin)-0.5 ) / 5. + 0.025; 
  }
};

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Function
FT sphere_function (const Point& p)
{ return CGAL::squared_distance(p, Point(CGAL::ORIGIN))-1; }

int main()
{
  // Domain (Warning: Sphere_3 constructor uses squared radius !)
  Mesh_domain domain(sphere_function,
                     K::Sphere_3(CGAL::ORIGIN, 2.));

  // Mesh criteria
  Spherical_sizing_field size;
  Mesh_criteria criteria(facet_angle=30, facet_size=0.1, facet_distance=0.025,
                         cell_radius_edge_ratio=2, cell_size=size);
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);

  return 0;
}

