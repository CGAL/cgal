#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/Real_timer.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT(Function)(const Point&);
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Remeshing
typedef CGAL::Triangulation_3<Tr::Geom_traits,
                              Tr::Triangulation_data_structure> T3_remeshing;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>     Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria  Cell_criteria;

// Sizing field
struct Spherical_sizing_field
{
  typedef ::FT FT;
  typedef Point Point_3;
  typedef Mesh_domain::Index Index;
  FT operator()(const Point_3& p, const int, const Index&) const
  {
    FT sq_d_to_origin = CGAL::squared_distance(p, Point(CGAL::ORIGIN));
    return CGAL::abs(CGAL::sqrt(sq_d_to_origin) - 0.5) / 5. + 0.025;
  }
};

// Function
FT sphere_function(const Point& p)
{
  return CGAL::squared_distance(p, Point(CGAL::ORIGIN)) - 1;
}

using namespace CGAL::parameters;

int main()
{
  CGAL::Real_timer timer;
  timer.start();

  Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain
                        (sphere_function, K::Sphere_3(CGAL::ORIGIN, K::FT(2)));

  // Mesh criteria
  Spherical_sizing_field size;
  Mesh_criteria criteria(facet_angle(30).facet_size(0.1).facet_distance(0.025).
                         cell_radius_edge_ratio(2).cell_size(size));

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude().no_perturb());

  timer.stop();
  std::cout << "Meshing done (" << timer.time()
            << " seconds)" << std::endl;

  //std::ofstream os_meshing("out_meshing.mesh");
  //CGAL::IO::write_MEDIT(os_meshing, c3t3.triangulation());
  //os_meshing.close();

  //Remeshing : extract triangulation
  timer.reset();
  timer.start();
  T3_remeshing t3 = CGAL::convert_to_triangulation_3(c3t3);

  //Remeshing : coarsen
  //double target_edge_length = 0.15;//for uniform
  CGAL::tetrahedral_isotropic_remeshing(t3, size,
      number_of_iterations(3));

  timer.stop();
  std::cout << "Remeshing done (" << timer.time()
            << " seconds)" << std::endl;

  //std::ofstream os_remeshing("out_remeshing.mesh");
  //CGAL::IO::write_MEDIT(os_remeshing, t3);
  //os_remeshing.close();

  return 0;
}
