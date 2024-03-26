#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_complex_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <cstdlib>
#include <cassert>

#include <fstream>
#include <CGAL/IO/File_medit.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_complex_mesh_domain_3<K> Mesh_domain;


typedef CGAL::Sequential_tag Concurrency_tag;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

namespace params = CGAL::parameters;

const char* const filenames[] = {
  "meshes/polyhedral_complex_of_spheres/Sphere1.off",
  "meshes/polyhedral_complex_of_spheres/Sphere2.off",
  "meshes/polyhedral_complex_of_spheres/Sphere3.off",
  "meshes/polyhedral_complex_of_spheres/Intersection12.off",
  "meshes/polyhedral_complex_of_spheres/Intersection13.off"
};

const std::pair<int, int> incident_subdomains[] = {
  std::make_pair(1, 0),
  std::make_pair(2, 0),
  std::make_pair(3, 0),
  std::make_pair(1, 2),
  std::make_pair(1, 3)
};

int main()
{
  const std::size_t nb_patches = sizeof(filenames) / sizeof(const char*);
  assert(sizeof(incident_subdomains) == nb_patches * sizeof(std::pair<int, int>));
  std::vector<Polyhedron> patches(nb_patches);
  for(std::size_t i = 0; i < nb_patches; ++i) {
    std::ifstream input(CGAL::data_file_path(filenames[i]));
    if(!(input >> patches[i])) {
      std::cerr << "Error reading " << filenames[i] << " as a polyhedron!/n";
      return EXIT_FAILURE;
    }
  }
  // Create domain
  Mesh_domain domain(patches.begin(), patches.end(),
                     incident_subdomains, incident_subdomains+nb_patches);

  // do not detect borders/features, otherwise the manifold criterion
  // will have no effect, and the test will be useless
  //domain.detect_borders();

  // Mesh criteria
  Mesh_criteria criteria(params::edge_size(1.0).edge_min_size(0.1)
    //.facet_distance(0.1)
    .facet_min_size(0.05)
    //.facet_topology(CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH)
    //.cell_radius_edge_ratio(3.)
    );

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
    params::no_perturb(), params::no_exude(),
    params::manifold());

  std::ofstream medit_file("out.mesh");
  CGAL::IO::write_MEDIT(medit_file, c3t3);
  medit_file.close();

  return EXIT_SUCCESS;
}
