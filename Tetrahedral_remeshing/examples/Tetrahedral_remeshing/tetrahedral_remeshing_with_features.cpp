#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/property_map.h>

#include "tetrahedral_remeshing_generate_input.h"

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Remeshing_triangulation = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;
using Point = Remeshing_triangulation::Point;
using Vertex_handle = Remeshing_triangulation::Vertex_handle;

using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
using Constraints_set = std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;


int main(int argc, char* argv[])
{
  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.02;
  const int nb_iter = (argc > 2) ? atoi(argv[2]) : 1;
  const int nbv = (argc > 3) ? atoi(argv[3]) : 500;

  Remeshing_triangulation t3;

  Constraints_set constraints;
  CGAL::Tetrahedral_remeshing::generate_input_cube(nbv, t3, constraints);
  make_constraints_from_cube_edges(t3, constraints);

  CGAL::tetrahedral_isotropic_remeshing(t3, target_edge_length,
    CGAL::parameters::edge_is_constrained_map(Constraints_pmap(constraints))
                     .number_of_iterations(nb_iter));

  return EXIT_SUCCESS;
}
