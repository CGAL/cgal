#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

#include <CGAL/tetrahedral_remeshing.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Image_3 Image;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;

typedef CGAL::Triangulation_3<typename Tr::Geom_traits,
                              typename Tr::Triangulation_data_structure> T3;

using namespace CGAL::parameters;

int main()
{
  const char* filename = "data/liver.inr.gz";

  CGAL::Image_3 image;
  if (!image.read(filename)) {
    std::cerr << "Error: Cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }
  Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image, 1e-9);

  // Mesh criteria
  Facet_criteria facet_criteria(25, 20, 2); // angle, size, approximation
  Cell_criteria cell_criteria(3, 20); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());

  std::cout << "Meshing done." << std::endl;

  //Remeshing : extract triangulation
  T3 t3 = CGAL::convert_to_triangulation_3(c3t3);

  //Remeshing : coarsen
  double target_edge_length = 15.;
  CGAL::tetrahedral_isotropic_remeshing(t3, target_edge_length,
    CGAL::parameters::number_of_iterations(2)
    .smooth_constrained_edges(true));

  std::cout << "Remeshing 1 done." << std::endl;

  //Remeshing : refine
  target_edge_length = 20.;
  CGAL::tetrahedral_isotropic_remeshing(t3, target_edge_length,
    CGAL::parameters::number_of_iterations(2).remesh_boundaries(false));

  std::cout << "Remeshing 2 done." << std::endl;

  return 0;
}
