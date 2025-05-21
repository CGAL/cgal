#define CGAL_NO_DEPRECATION_WARNINGS

#include "test_dependencies.h"
#include <CGAL/config.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>

#include <CGAL/lloyd_optimize_mesh_2.h>
#include <CGAL/Mesh_optimization_return_code.h>

#include "test_utilities.h"

#include <CGAL/IO/File_poly.h>
#include <iostream>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>  Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>    Fb;

typedef CGAL::Triangulation_data_structure_2<Vb, Fb>  TDS;
typedef CGAL::Exact_predicates_tag                    Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::size_type size_type;
typedef CDT::Point Point;

using namespace CGAL::parameters;

struct Lloyd_tester
{
  void operator()(CDT& cdt) const
  {
    std::vector<Point> seeds;
    seeds.reserve(32);

    std::cerr << "Reading fish-and-rectangle.poly...";
    std::ifstream poly_file("fish-and-rectangle.poly");
    CGAL::IO::read_triangle_poly_file(cdt, poly_file, std::back_inserter(seeds));
    assert( cdt.is_valid() );

    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\nNumber of seeds: " << seeds.size() << "\n\n";

    std::cerr << "Meshing the triangulation with size 0.1...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 CGAL::parameters::seeds(seeds).criteria(Criteria(0.125, 0.1)));
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() << "\n\n";
    assert( cdt.is_valid() );
    assert( 580 <= cdt.number_of_vertices() &&
                    cdt.number_of_vertices() <= 640 );


    const size_type number_of_constraints = number_of_constrained_edges(cdt);
    const size_type number_of_vertices1 = cdt.number_of_vertices();

    CGAL::Mesh_optimization_return_code rc
      = CGAL::lloyd_optimize_mesh_2(cdt,
              CGAL::parameters::number_of_iterations(10).
              convergence(0.001).
              freeze_bound(0.001).
              seeds(seeds));
    const size_type number_of_vertices2 = cdt.number_of_vertices();
    std::cerr << " done (return code = "<< rc <<").\n";
    std::cerr << "Number of vertices: " << number_of_vertices2 << "\n\n";

    assert( cdt.is_valid() );
    assert( number_of_vertices1 == number_of_vertices2 );
    assert( number_of_constraints == number_of_constrained_edges(cdt));
  }
};

struct Lloyd_tester_original_BP_API
{
  void operator()(CDT& cdt) const
  {
    std::vector<Point> seeds;
    seeds.reserve(32);

    std::cerr << "Reading fish-and-rectangle.poly...";
    std::ifstream poly_file("fish-and-rectangle.poly");
    CGAL::IO::read_triangle_poly_file(cdt, poly_file, std::back_inserter(seeds));
    assert( cdt.is_valid() );

    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\nNumber of seeds: " << seeds.size() << "\n\n";

    std::cerr << "Meshing the triangulation with size 0.1...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 CGAL::parameters::seeds(seeds).
                                 criteria(Criteria(0.125, 0.1)));
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() << "\n\n";
    assert( cdt.is_valid() );
    assert( 580 <= cdt.number_of_vertices() &&
                    cdt.number_of_vertices() <= 640 );


    const size_type number_of_constraints = number_of_constrained_edges(cdt);
    const size_type number_of_vertices1 = cdt.number_of_vertices();

    CGAL::Mesh_optimization_return_code rc
      = CGAL::lloyd_optimize_mesh_2(cdt,
              max_iteration_number = 10,
              convergence = 0.001,
              freeze_bound = 0.001,
              seeds_begin = seeds.begin(),
              seeds_end = seeds.end());
    const size_type number_of_vertices2 = cdt.number_of_vertices();
    std::cerr << " done (return code = "<< rc <<").\n";
    std::cerr << "Number of vertices: " << number_of_vertices2 << "\n\n";

    assert( cdt.is_valid() );
    assert( number_of_vertices1 == number_of_vertices2 );
    assert( number_of_constraints == number_of_constrained_edges(cdt));
  }
};

int main()
{
  std::cerr << "TESTING lloyd_optimize_mesh_2 with Epick...\n\n";
  CDT cdt;
  Lloyd_tester tester;
  tester(cdt);

  std::cerr << "TESTING lloyd_optimize_mesh_2 with Epick (original Boost Parameter API)...\n\n";
  cdt = CDT();
  Lloyd_tester_original_BP_API tester_bis;
  tester_bis(cdt);

  return 0;
}
