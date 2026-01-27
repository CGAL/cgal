#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_triangulation_3.h>

#include <CGAL/Image_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/IO/File_medit.h>

#include <string>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Image_3 Image;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

typedef Tr::Geom_traits Gt;
typedef CGAL::Triangulation_3<Gt, typename Tr::Triangulation_data_structure> T3;
typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<Gt> Remeshing_triangulation;

using namespace CGAL::parameters;

int main()
{
  const std::string filename = CGAL::data_file_path("images/liver.inr.gz");

  CGAL::Image_3 image;
  if(!image.read(filename)) {
    std::cerr << "Error: Cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }
  Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image, relative_error_bound = 1e-9);

  // Mesh criteria
  Facet_criteria facet_criteria(25, 20, 2); // angle, size, approximation
  Cell_criteria cell_criteria(3, 20);       // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());

  std::cout << "Meshing done." << std::endl;

  // Write
  std::ofstream os("after_meshing_io.mesh");
  CGAL::IO::write_MEDIT(os, c3t3);
  os.close();

  // Read
  Remeshing_triangulation tr;
  std::ifstream is("after_meshing_io.mesh");
  CGAL::IO::read_MEDIT(is, tr);

  // Remesh
  double target_edge_length = 10.;
  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length);

  std::cout << "Remeshing done." << std::endl;

  return 0;
}
