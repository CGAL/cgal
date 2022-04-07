// 154 515 565
#include <CGAL/config.h>
#include "test_dependencies.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#if CGAL_USE_CORE || CGAL_USE_LEDA
#  include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#endif
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_no_edge_refinement_2.h>

#include <CGAL/IO/File_poly.h>

#include <fstream>
#include <iostream>
#include <cassert>

#include "test_utilities.h"

template <typename K, typename CDT>
struct Tester2;

template <typename K>
struct Tester {
  void operator()() const {

    typedef CGAL::Triangulation_vertex_base_2<K> Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;

    typedef CGAL::Constrained_triangulation_plus_2<CDT> CDTplus;

    std::cerr << "  ### Testing With CDT...\n";
    Tester2<K, CDT> tester;
    tester();
    std::cerr << "  ### Testing With CDT_plus_2...\n";
    Tester2<K, CDTplus> tester_plus;
    tester_plus();
  }
};

template <typename K, typename CDT>
struct Tester2 {
  void operator()() const {
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
    typedef typename CDT::size_type size_type;

    typedef typename CDT::Point Point;

    CDT cdt;

    std::vector<Point> seeds;
    seeds.reserve(32);

    std::cerr << "Reading fish-and-rectangle.poly...";
    std::ifstream poly_file("fish-and-rectangle.poly");
    CGAL::IO::read_triangle_poly_file(cdt, poly_file, std::back_inserter(seeds));
    assert(cdt.is_valid());
    const size_type inititial_number_of_vertices = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\nNumber of seeds: " << seeds.size() << "\n\n";

    std::cerr << "Saving the triangulation...\n\n";
    CDT cdt2 = cdt;

    std::cerr << "1/ First tests:\n\n";

    std::cerr << "Meshing the triangulation with size 0...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 seeds.begin(), seeds.end(),
                                 Criteria());
    const size_type number_of_vertices0 = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() << "\n\n";
    assert( 64 <= cdt.number_of_vertices() &&
                    cdt.number_of_vertices() <= 72 );
    assert( seeds.size() == 3 );

    std::cerr << "Meshing the triangulation with size 0.2...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 seeds.begin(), seeds.end(),
                                 Criteria(0.125, 0.2));

    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() << "\n\n";
    assert( 190 <= cdt.number_of_vertices() &&
                    cdt.number_of_vertices() <= 210 );

    std::cerr << "Meshing the triangulation with size 0.1...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 seeds.begin(), seeds.end(),
                                 Criteria(0.125, 0.1));
    const size_type number_of_vertices1 = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() << "\n\n";
    assert( 580 <= cdt.number_of_vertices() &&
                    cdt.number_of_vertices() <= 640 );

    cdt = cdt2;
    std::cerr << "Triangulation restored.\n";
    std::cerr << "Number of vertices: " << cdt.number_of_vertices() << "\n\n";

    std::cerr << "Meshing the triangulation with Delaunay_mesh_criteria_2<CDT>()...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 seeds.begin(), seeds.end(),
                                 CGAL::Delaunay_mesh_criteria_2<CDT>());
    const size_type number_of_vertices0bis = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() << "\n\n";

    assert( number_of_vertices0 == number_of_vertices0bis );

    cdt = cdt2;
    std::cerr << "Triangulation restored.\n";
    std::cerr << "Number of vertices: " << cdt.number_of_vertices() << "\n\n";

    std::cerr << "2/ Comparaison between refine_Delaunay_mesh_2() and other"
              << " possibilities:\n\n";

    std::cerr << "Meshing the triangulation with size 0.1, with "
              << "refine_Delaunay_mesh_2()...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 seeds.begin(), seeds.end(),
                                 Criteria(0.125, 0.1));
    const size_type number_of_vertices1bis = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\n\n";

    assert( number_of_vertices1bis <= number_of_vertices1 );

    cdt = cdt2;

    std::cerr << "Meshing the triangulation with size 0.1, with "
              << "mesher.refine_mesh()...";
    {
      CGAL::Delaunay_mesher_2<CDT, Criteria> mesher(cdt, Criteria(0.125, 0.1));
      mesher.set_seeds(seeds.begin(), seeds.end());
      mesher.refine_mesh();
    }
    const size_type number_of_vertices2 = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\n\n";

    assert(cdt.is_valid());
    assert( number_of_vertices2 == number_of_vertices1bis );

    cdt = cdt2;

    std::cerr << "Meshing the triangulation with size 0.1, with\n"
              << "a loop of mesher.try_one_step_refine_mesh()...";
    size_type step = 0;
    {
      CGAL::Delaunay_mesher_2<CDT, Criteria> mesher(cdt, Criteria(0.125, 0.1));
      mesher.set_seeds(seeds.begin(), seeds.end());
      mesher.init();
      while(mesher.try_one_step_refine_mesh())
        ++step;
    }
    const size_type number_of_vertices3 = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\nNumber of steps: " << step << "\n\n";

    assert( step + inititial_number_of_vertices >= number_of_vertices3 );
    assert( number_of_vertices3 == number_of_vertices2 );

    cdt = cdt2;

    std::cerr << "Meshing the triangulation with size 0.1, with\n"
              << "a loop of mesher.step_by_step_refine_mesh()...";
    step = 0;
    {
      CGAL::Delaunay_mesher_2<CDT, Criteria> mesher(cdt, Criteria(0.125, 0.1));
      mesher.set_seeds(seeds.begin(), seeds.end());
      mesher.init();
      while(mesher.step_by_step_refine_mesh())
        ++step;
    }
    const size_type number_of_vertices4 = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\nNumber of steps: " << step << "\n\n";

    assert( number_of_vertices4 == number_of_vertices2 );
    assert( number_of_vertices4 == step + inititial_number_of_vertices );

    std::cerr << "Test the undocumented function:"
              << "  refine_Delaunay_mesh_2_without_edge_refinement\n"
              << "with size 0.1...";
    cdt = cdt2;
    CGAL::refine_Delaunay_mesh_2_without_edge_refinement(cdt, Criteria(0.125, 0.1));
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\n";
    assert(cdt.number_of_vertices() == 36);
  }
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K_e_i;
#if CGAL_USE_CORE || CGAL_USE_LEDA
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt K_e_e;
#endif

int main()
{
  std::cerr << "TESTING WITH Exact_predicates_inexact_constructions_kernel...\n\n";
  Tester<K_e_i> tester1;
  tester1();
#if CGAL_USE_CORE || CGAL_USE_LEDA
  // std::cerr << "\n\nTESTING WITH Exact_predicates_exact_constructions_kernel_with_sqrt...\n\n";
  // Tester<K_e_e> tester2;
  // tester2();
#endif
}
