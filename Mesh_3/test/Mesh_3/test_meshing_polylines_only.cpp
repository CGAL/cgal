// #define CGAL_MESH_3_PROTECTION_HIGH_VERBOSITY 1
// #define CGAL_MESH_3_PROTECTION_DEBUG 255
#define CGAL_PROFILE 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Intersections_3/Sphere_3_Sphere_3.h>
#include <CGAL/IO/File_binary_mesh_3.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/iterator.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_3/Protect_edges_sizing_field.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Random.h>
#include <CGAL/Sizing_field_with_aabb_tree.h>
#include <CGAL/Sphere_3.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

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
  std::cout << "\tSeed is\t" << CGAL::get_default_random().get_seed() << std::endl;
  std::cerr.precision(17);
  std::cout.precision(17);

  std::string input_path = (argc > 1) ? argv[1] : CGAL::data_file_path("polylines_3/couplingdown-polylines.txt");
  std::string output_path = (argc > 2) ? argv[2] : "";

  typedef K::Point_3 Point;

  Polyhedron p;
  bool file_is_a_polyhedron = CGAL::IO::read_polygon_mesh(input_path, p);
  if(file_is_a_polyhedron) {
    std::cout << "Read polyhedron from " << input_path << std::endl;
  } else {
    // Create domain with a fake polyhedron
    p.clear();
    p.make_tetrahedron(Point(0, 0, 0),
                      Point(1, 0, 0),
                      Point(0, 1, 0),
                      Point(0, 0, 1));
  }

  Mesh_domain domain(p, &CGAL::get_default_random());

  if(file_is_a_polyhedron) {
    std::cout << "Extracting features from polyhedron...\n";
    domain.detect_features();
  } else {
    typedef std::vector<K::Point_3> Polyline;
    typedef std::vector<Polyline> Polylines;

    Polylines polylines;
    std::ifstream in(input_path);
    while(!in.eof()) {
      Polyline polyline;
      std::size_t n;
      if(!(in >> n)) {
        if(in.eof()) continue;
        else return 1;
      }
      // std::cerr << "Reading polyline #" << polylines.size()
      //           << " with " << n << " vertices\n";
      polyline.reserve(n);
      while( n > 0 ) {
        K::Point_3 p;
        if(!(in >> p)) return 1;
        polyline.push_back(p);
        --n;
      }
      polylines.push_back(polyline);
    }
    std::cerr << "Read " << polylines.size() << " polylines from " << input_path << std::endl;
    domain.add_features(polylines.begin(), polylines.end());
  }

  std::size_t n = 0;
  auto out = boost::make_function_output_iterator([&n](auto&&...) { ++n; });
  domain.get_corners(out);
  std::cout << "Number of corners: " << n << std::endl;
  n = 0;
  domain.get_curves(out);
  std::cout << "Number of curves: " << n << std::endl;

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.1);
  typedef Mesh_criteria::Edge_criteria Edge_criteria;

  C3t3 c3t3;
  typedef CGAL::Mesh_3::internal::Edge_criteria_sizing_field_wrapper<Edge_criteria> Sizing_field;
  CGAL::Mesh_3::Protect_edges_sizing_field<C3t3, Mesh_domain, Sizing_field>
    protect_edges(c3t3, domain, Sizing_field(criteria.edge_criteria_object()), 0.01);
  protect_edges(true);

  auto& tr = c3t3.triangulation();
  std::cout << "Number of vertices in c3t3: "
            << tr.number_of_vertices() << std::endl;
  std::cout << "Number of corners in c3t3: "
            << c3t3.number_of_corners() << std::endl;
  std::cout << "Number of edges in c3t3: "
            << c3t3.number_of_edges() << std::endl;

  if(!output_path.empty()) {
    std::ofstream out(output_path);
    out.precision(17);
    for(auto v: tr.finite_vertex_handles()) {
      out << v->point().point() << ' ' << CGAL::sqrt(v->point().weight()) << '\n';
    }
  }

  assert(c3t3.is_valid());

  // the following std::transform_reduce is like std::all_of but without the short-circuiting
  bool ok = std::transform_reduce(tr.finite_edges_begin(), tr.finite_edges_end(), true, std::logical_and<>{},
    [&](auto e) {
      auto [v, w] = tr.vertices(e);
      assert(v->in_dimension() >= 0);
      assert(v->in_dimension() <= 1);
      assert(w->in_dimension() >= 0);
      assert(w->in_dimension() <= 1);
      if(v->is_special() || w->is_special()) {
        return true;
      }
      auto curve_index_of_e = c3t3.curve_index(e);
      bool e_is_in_complex = (curve_index_of_e != C3t3::Curve_index());
      const auto& p1 = tr.point(v);
      const auto& p2 = tr.point(w);
      bool balls_intersect =
          do_intersect(CGAL::Sphere_3<K>(p1.point(), p1.weight()), CGAL::Sphere_3<K>(p2.point(), p2.weight()));
      bool e_ok = e_is_in_complex == balls_intersect;
      if(!e_ok) {
        std::cerr << "Error: two " << (e_is_in_complex ? "adjacent" : "non-adjacent") << " vertices have "
                  << (balls_intersect ? "intersecting" : "non-intersecting") << " protection balls\n";
        std::cerr << "v: " << p1 << " with radius " << CGAL::sqrt(p1.weight()) << "\n";
        std::cerr << "w: " << p2 << " with radius " << CGAL::sqrt(p2.weight()) << "\n";
      }
      if(e_is_in_complex) {
        for(auto v : tr.vertices(e)) {
          if(v->in_dimension() == 1) {
            auto v_curve_index = domain.curve_index(v->index());
            if(v_curve_index != curve_index_of_e) {
              std::cerr << "Error: vertex in dimension 1 is incident to an edge in the complex, but not on the same curve.\n";
              std::cerr << "v: " << p1 << " with radius " << CGAL::sqrt(p1.weight()) << "\n";
              std::cerr << "e curve index: " << curve_index_of_e << "\n";
              std::cerr << "v curve index: " << v_curve_index << "\n";
              e_ok = false;
            }
          }
        }
      }
      return e_ok;
    }
  );
  if(!ok) {
    return EXIT_FAILURE;
  }
}
