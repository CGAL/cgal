#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef CGAL::Surface_mesh<K::Point_3>                          Mesh;

typedef boost::graph_traits<Mesh>::edge_descriptor              edge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const char* filename = argc > 1 ? argv[1] : "data/anchor_dense.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid .off file." << std::endl;
    return EXIT_FAILURE;
  }

  // Constrain edges with a dihedral angle over 60°
  typedef boost::property_map<Mesh, CGAL::edge_is_feature_t>::type EIFMap;
  EIFMap eif = get(CGAL::edge_is_feature, mesh);
  PMP::detect_sharp_edges(mesh, 60, eif);

  int sharp_counter = 0;
  for(edge_descriptor e : edges(mesh))
    if(get(eif, e))
      ++sharp_counter;

  std::cout << sharp_counter << " sharp edges" << std::endl;

  const unsigned int nb_iterations = (argc > 2) ? std::atoi(argv[2]) : 10;

  std::cout << "Smoothing mesh... (" << nb_iterations << " iterations)" << std::endl;

  // Smooth with both angle and area criteria + Delaunay flips
  PMP::smooth_mesh(mesh, PMP::parameters::number_of_iterations(nb_iterations)
                                         .use_safety_constraints(false) // authorize all moves
                                         .edge_is_constrained_map(eif));

  std::ofstream output("mesh_smoothed.off");
  output.precision(17);
  output << mesh;

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
