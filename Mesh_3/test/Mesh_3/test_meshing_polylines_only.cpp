//#define CGAL_MESH_3_PROTECTION_DEBUG 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/IO/File_binary_mesh_3.h>

#include <cassert>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;


#include <fstream>

int main(int argc, char** argv)
{
  if(argc != 2) {
    std::cerr << "This test needs a filename as argument.\n";
    return 1;
  }
  typedef K::Point_3 Point;

  // Create domain
  Polyhedron p;
  p.make_tetrahedron(Point(0, 0, 0),
                     Point(1, 0, 0),
                     Point(0, 1, 0),
                     Point(0, 0, 1));

    std::cout << "\tSeed is\t"
      << CGAL::get_default_random().get_seed() << std::endl;
  Mesh_domain domain(p, &CGAL::get_default_random());

  typedef std::vector<K::Point_3> Polyline;
  typedef std::vector<Polyline> Polylines;

  C3t3 c3t3;

  Polylines polylines;
  std::ifstream in(argv[1]);
  while(!in.eof()) {
    Polyline polyline;
    std::size_t n;
    if(!(in >> n)) {
      if(in.eof()) continue;
      else return 1;
    }
    std::cerr << "Reading polyline #" << polylines.size()
              << " with " << n << " vertices\n";
    polyline.reserve(n);
    while( n > 0 ) {
      K::Point_3 p;
      if(!(in >> p)) return 1;
      polyline.push_back(p);
      --n;
    }
    polylines.push_back(polyline);
  }
  std::cerr << "Number of polylines: " << polylines.size() << std::endl;

  domain.add_features(polylines.begin(), polylines.end());

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.1);
  typedef Mesh_criteria::Edge_criteria Edge_criteria;
  typedef CGAL::Mesh_3::internal::Edge_criteria_sizing_field_wrapper<Edge_criteria> Sizing_field;
  CGAL::Mesh_3::Protect_edges_sizing_field<C3t3, Mesh_domain, Sizing_field>
    protect_edges(c3t3, domain, Sizing_field(criteria.edge_criteria_object()), 0.01);
  protect_edges(true);

  // CGAL::Mesh_3::internal::init_c3t3_with_features(c3t3, domain, criteria);

  // Output
  std::ofstream medit_file("out-mesh-polylines.mesh");
  c3t3.output_to_medit(medit_file);
  std::ofstream binary_file("out-mesh-polylines.binary.cgal", std::ios::binary|std::ios::out);
  CGAL::IO::save_binary_file(binary_file, c3t3);
  std::cout << "Number of vertices in c3t3: "
            << c3t3.triangulation().number_of_vertices() << std::endl;
  assert(c3t3.triangulation().number_of_vertices() > 900);
  assert(c3t3.triangulation().number_of_vertices() < 1100);
}
