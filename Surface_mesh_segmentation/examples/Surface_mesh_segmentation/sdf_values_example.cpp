#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>

#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main()
{
  // create and read Polyhedron
  Polyhedron mesh;
  std::ifstream input("data/cactus.off");
  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Not a valid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // create a property-map
  typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
  Facet_double_map internal_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(internal_map);

  // compute SDF values
  std::pair<double, double> min_max_sdf = CGAL::sdf_values(mesh, sdf_property_map);


  // It is possible to compute the raw SDF values and post-process them using
  // the following lines:
  // const std::size_t number_of_rays = 25;  // cast 25 rays per facet
  // const double cone_angle = 2.0 / 3.0 * CGAL_PI; // set cone opening-angle
  // CGAL::sdf_values(mesh, sdf_property_map, cone_angle, number_of_rays, false);
  // std::pair<double, double> min_max_sdf =
  //  CGAL::sdf_values_postprocessing(mesh, sdf_property_map);

  // print minimum & maximum SDF values
  std::cout << "minimum SDF: " << min_max_sdf.first
            << " maximum SDF: " << min_max_sdf.second << std::endl;

  // print SDF values
  for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin();
      facet_it != mesh.facets_end(); ++facet_it) {
      std::cout << sdf_property_map[facet_it] << " ";
  }
  std::cout << std::endl;
}
