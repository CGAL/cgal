#define CGAL_MESH_3_USE_EXPERIMENTAL_COMPACT_MESH_CELL_BASE_3
#define CGAL_MESH_3_NO_CIRCUMCENTER_CACHE

#include "test_meshing_utilities.h"
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/File_tetgen.h>
#include <CGAL/IO/File_binary_mesh_3.h>
#include <CGAL/use.h>

#include <fstream>
#include <sstream>

static constexpr bool verbose =
#ifdef CGAL_MESH_3_VERBOSE
  true;
#else
  false;
#endif


std::stringstream cerr_output;

template <typename K, typename Concurrency_tag = CGAL::Sequential_tag>
struct Polyhedron_with_features_tester
  : public Tester<K>
{
  std::streambuf* old_cerr_buf;
  std::string tag_name;

  Polyhedron_with_features_tester(std::string tag_name)
    : old_cerr_buf(std::cerr.rdbuf(verbose ? cerr_output.rdbuf() : std::cerr.rdbuf()))
    , tag_name(tag_name)
  {
  }
  ~Polyhedron_with_features_tester()
  {
    std::cerr.rdbuf(old_cerr_buf);
    if(verbose) {
      const auto output = cerr_output.str();
      const auto str_size= output.size();
      const auto str_begin = output.data();
      const auto str_end = str_begin + str_size;
      constexpr auto nb = static_cast<std::remove_cv_t<decltype(str_size)>>(10000);
      auto pos1 = str_begin + (std::min)(nb, str_size);
      assert(pos1 <= str_end);
      const auto pos2 = str_end - (std::min)(nb, str_size);
      assert(pos2 >= str_begin);
      if(pos2 <= pos1) {
        pos1 = str_end;
      }
      std::cerr << "NOW THE FIRST AND LAST 10k CHARACTERS OF THE CERR OUTPUT:\n";
      std::cerr << "-----\n" << std::string(str_begin, pos1) << "\n-----\n";
      if(pos1 != str_end) {
        std::cerr << "[...]\n-----\n" << std::string(pos2, str_end) << "\n-----\n";
      }
      const auto file_name = std::string("test_meshing_verbose-") + tag_name + ".txt";
      std::ofstream file_output(file_name);
      file_output << output;
      std::cerr << "Full log is output to " << file_name << "\n";
    }
  }
  void operator()() const
  {
    typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> GT;
    typedef typename CGAL::Mesh_polyhedron_3<GT, short>::type Polyhedron;
    typedef CGAL::Polyhedral_mesh_domain_with_features_3<GT,
                                                         Polyhedron,
                                                         CGAL::Default,
                                                         short> Mesh_domain;

    typedef typename CGAL::Mesh_triangulation_3<
      Mesh_domain,
      typename CGAL::Kernel_traits<Mesh_domain>::Kernel,
      Concurrency_tag>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3 <
      Tr,
      typename Mesh_domain::Corner_index,
      typename Mesh_domain::Curve_index > C3t3;

    typedef CGAL::Mesh_criteria_3<C3t3> Mesh_criteria;
    typedef typename Mesh_criteria::Edge_criteria Edge_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    std::cout << "\tSeed is\t" << CGAL::get_default_random().get_seed() << std::endl;

    std::ifstream input(CGAL::data_file_path("meshes/cube.off"));
    Polyhedron polyhedron;
    input >> polyhedron;

    Mesh_domain domain(polyhedron, &CGAL::get_default_random());
    domain.detect_features();

    // Set mesh criteria
#ifndef CGAL_MESH_3_VERBOSE
    Edge_criteria edge_criteria(0.2);
    Facet_criteria facet_criteria(30, 0.2, 0.02);
    Cell_criteria cell_criteria(3, 0.2);
#else // a different set of criteria, for the test of CGAL_MESH_3_VERBOSE
    Edge_criteria edge_criteria(0.3);
    Facet_criteria facet_criteria(30, 0.3, 0.03);
    Cell_criteria cell_criteria(3, 0.4);
#endif
    Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

    namespace params = CGAL::parameters;

    // Mesh generation
    std::cout << "sizeof(Vertex) " << sizeof(typename C3t3::Triangulation::Vertex) << std::endl;
    std::cout << "sizeof(Cell) " << sizeof(typename C3t3::Triangulation::Cell) << std::endl;

    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        params::odt(params::convergence(0.03).freeze_bound(0.02).time_limit(30)),
                                        params::lloyd(params::max_iteration_number(10)),
                                        params::perturb(params::sliver_bound(10).time_limit(30)),
                                        params::exude(params::sliver_bound(10).time_limit(0)));

    CGAL::remove_far_points_in_mesh_3(c3t3);
  }
};

int main()
{
  {
    std::cerr << "Mesh generation from a polyhedron with edges:\n";
    Polyhedron_with_features_tester<K_e_i> test_epic("sequential");
    test_epic();
  }
#ifdef CGAL_LINKED_WITH_TBB
  {
    std::cerr << "\n\nParallel mesh generation from a polyhedron with edges:\n";
    Polyhedron_with_features_tester<K_e_i, CGAL::Parallel_tag> test_epic_p("parallel");
    test_epic_p();
  }
#endif

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
