#define INDEX_STORAGE  1
//#define CGAL_MESH_3_VERBOSE
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Memory_sizer.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

typedef CGAL::Memory_sizer Memory_sizer;

namespace params = CGAL::parameters;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("meshes/elephant.off");
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input(fname);
  input >> polyhedron;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  input.close();

  if (!CGAL::is_triangle_mesh(polyhedron)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // Create domain
  Mesh_domain domain(polyhedron);

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(params::facet_angle(25).facet_size(0.15).facet_distance(0.008).
                                 cell_radius_edge_ratio(3));

  Memory_sizer memory_sizer;
  auto res_mem_at_start = memory_sizer.resident_size();
  std::cout << "Memory usage before CGAL::make_mesh_3:\n" << memory_sizer.virtual_size() << " bytes (virtual), "
            << res_mem_at_start << " bytes (resident)" << std::endl;

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, params::no_perturb().no_exude());

  auto res_mem = memory_sizer.resident_size();
  std::cout << "Memory usage after construction of the triangulation:\n" << memory_sizer.virtual_size() << " bytes (virtual), "
            << res_mem << " bytes (resident)" << std::endl;
  std::cout << "Diff in resident memory: "
            << res_mem - res_mem_at_start << " bytes" << std::endl;

  std::cout << "cell size  : " << c3t3.triangulation().infinite_vertex()->cell()->size() << std::endl;
  std::cout << "vertex size: " << c3t3.triangulation().infinite_vertex()->size() << std::endl;

  // Output
  std::ofstream medit_file("out_1.mesh");
  CGAL::IO::write_MEDIT(medit_file, c3t3);
  medit_file.close();

  // Set tetrahedron size (keep cell_radius_edge_ratio), ignore facets
  Mesh_criteria new_criteria(params::cell_radius_edge_ratio(3).cell_size(0.03));

  // Mesh refinement
  CGAL::refine_mesh_3(c3t3, domain, new_criteria);

  // Output
  medit_file.open("out_2.mesh");
  CGAL::IO::write_MEDIT(medit_file, c3t3);
  medit_file.close();

  return EXIT_SUCCESS;
}
