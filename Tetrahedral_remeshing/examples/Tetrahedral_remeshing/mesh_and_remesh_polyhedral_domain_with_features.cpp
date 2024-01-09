#define CGAL_MESH_3_VERBOSE 1
#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/IO/File_binary_mesh_3.h>

#include <boost/container/flat_set.hpp>

int nb_surface_flip_candidates = 0;
int nb_surface_flip_done = 0;
int nb_surface_44_configs = 0;
int nb_surface_nm_configs = 0;
int nb_surface_44_flips_done = 0;
int nb_surface_nm_flips_done = 0;

std::ostringstream oss_flip("flipped_surface_edges.polylines.txt");


// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron> Mesh_domain;


#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char* argv[])
{
  //"data/Shape1/shape1.off" 40. "data/Shape1/mesh_40.binary.cgal" "data/Shape1/remesh_40.binary.cgal"
  const char* fname = (argc > 1) ? argv[1] : "data/tensileASCII.off";
  const double target = (argc > 2) ? atof(argv[2]) : 200.;
  const int nb_iter = (argc > 3) ? atoi(argv[3]) : 10;
  const char* fname_mesh3 = (argc > 4) ? argv[4] : "data/Tensile/mesh_200.binary.cgal";
  const char* fname_remesh = (argc > 5) ? argv[5] : "data/Tensile/remesh_200.binary.cgal";

  std::ifstream input(fname);
  Polyhedron polyhedron;

  std::string filename(fname);
  //  if (filename.substr(filename.find_last_of(".")).compare(".off") == 0)
  input >> polyhedron;
  if (input.fail()) {
    std::cerr << "Error: Cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  if (!CGAL::is_triangle_mesh(polyhedron)) {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // Create domain
  Mesh_domain domain(polyhedron);

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  double CellSize = 31.3082;
  int mesh_factor = 50;
  double size = (CellSize * mesh_factor) / 100;
  //const double size = 40.;
  Mesh_criteria criteria(edge_size = size,
    facet_angle = 25,
    facet_size = size,
    cell_radius_edge_ratio = 2,
    cell_size = size);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);// , no_perturb(), no_exude());

  // Output
//  std::ofstream os(fname_mesh3, std::ios_base::out | std::ios_base::binary);
//  CGAL::Mesh_3::save_binary_file(os, c3t3);
  CGAL::dump_c3t3(c3t3, "out_after_meshing");

  const double min_dihedral_angle = CGAL::Tetrahedral_remeshing::min_dihedral_angle(c3t3);
  std::cout << "Min dihedral angle after meshing = " << min_dihedral_angle << std::endl;

  // Remeshing
  CGAL::tetrahedral_isotropic_remeshing(c3t3, size,
    CGAL::parameters::number_of_iterations(nb_iter));
  //    .smooth_constrained_edges(true));

    // Output
  //  std::ofstream osr(fname_remesh, std::ios_base::out | std::ios_base::binary);
  //  CGAL::Mesh_3::save_binary_file(osr, c3t3);
  CGAL::dump_c3t3(c3t3, "out_after_remeshing");

  std::cout << "Meshing and remeshing done." << std::endl;

  std::cout << "Surface flip attempts : " << nb_surface_flip_candidates << std::endl;
  std::cout << "Surface flip done     : " << nb_surface_flip_done << std::endl;
  std::cout << "Surface 4/4 flips candidates proportion : "
    << 100.*nb_surface_44_configs / (double)(nb_surface_nm_configs + nb_surface_44_configs)
    << " percent" << std::endl;

  std::cout << "nb surface flips 4/4 done = " << nb_surface_44_flips_done << std::endl;
  std::cout << "nb surface flips n/m done = " << nb_surface_nm_flips_done << std::endl;


  return EXIT_SUCCESS;
}
