#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/tetrahedral_remeshing.h>


#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_if_available_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// Triangulation for Remeshing
using T3 = CGAL::Triangulation_3<Tr::Geom_traits,
                                 Tr::Triangulation_data_structure>;

namespace params = CGAL::parameters;

// Sizing field
struct Spherical_sizing_field
{
  typedef K::FT FT;
  typedef K::Point_3 Point_3;
  typedef Mesh_domain::Index Index;

  FT operator()(const Point_3& p, const int, const Index&) const
  {
    FT sq_d_to_origin = CGAL::squared_distance(p, Point(CGAL::ORIGIN));
    return CGAL::abs(CGAL::sqrt(sq_d_to_origin) - 0.5) / 5. + 0.025;
  }
};

// Function
FT sphere_function (const Point& p)
{
  return CGAL::squared_distance(p, Point(CGAL::ORIGIN)) - 1;
}

int main()
{
  /// [Domain creation] (Warning: Sphere_3 constructor uses squared radius !)
  Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain
  (sphere_function, K::Sphere_3(CGAL::ORIGIN, K::FT(2))
  );
  /// [Domain creation]

  // Mesh criteria
  Spherical_sizing_field size;
  Mesh_criteria criteria(params::facet_angle(30).facet_size(0.1).facet_distance(0.025).
                         cell_radius_edge_ratio(2).cell_size(size));

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, params::no_exude().no_perturb());

  T3 tr = CGAL::convert_to_triangulation_3(std::move(c3t3));
  //note we use the move semantic, with std::move(c3t3),
  //  to avoid a copy of the triangulation by the function
  //  `CGAL::convert_to_triangulation_3()`
  //  After the call to this function, c3t3 is an empty and valid C3t3.
  //It is possible to use :  CGAL::convert_to_triangulation_3(c3t3),
  //  Then the triangulation is copied and duplicated, and c3t3 remains as is.

  std::cout << "Remeshing...";
  std::cout.flush();

  CGAL::tetrahedral_isotropic_remeshing(tr, size);

  std::cout << "\rRemeshing done." << std::endl;

  return EXIT_SUCCESS;
}
