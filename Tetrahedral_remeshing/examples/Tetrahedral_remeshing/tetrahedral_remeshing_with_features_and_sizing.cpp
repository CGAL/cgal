#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#define CGAL_TETRAHEDRAL_REMESHING_DEBUG
#define CGAL_DUMP_REMESHING_STEPS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/property_map.h>

#include <unordered_set>
#include <iostream>
#include <utility>
#include <cassert>

#include "tetrahedral_remeshing_generate_input.h"

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Remeshing_triangulation = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;
using Point = Remeshing_triangulation::Point;
using Vertex_handle = Remeshing_triangulation::Vertex_handle;

using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
using Constraints_set = std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;


struct Distance_from_corner_sizing_field
{
  typedef K::FT FT;
  typedef Point Point_3;

  const Point_3 corner = Point_3{-2., -2., -2}; //lower corner of the cube

  template<typename Index>
  FT operator()(const Point_3& p, const int, const Index&) const
  {
    const FT d_to_origin = CGAL::approximate_sqrt(CGAL::squared_distance(p, corner));
    return 0.02 + 0.1 * d_to_origin;
  }
};


int main(int argc, char* argv[])
{
  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.02;
  const int nb_iter = (argc > 2) ? atoi(argv[2]) : 1;
  const int nbv = (argc > 3) ? atoi(argv[3]) : 50;

  Remeshing_triangulation t3;

  Constraints_set constraints;
  CGAL::Tetrahedral_remeshing::generate_input_cube(nbv, t3, constraints);
  set_subdomain_for_finite_cells(t3, 1);
  set_surface_patch_for_convex_hull_facets(t3, 2);

  std::ofstream ofs("out0.mesh");
  CGAL::IO::write_MEDIT(ofs, t3, CGAL::parameters::rebind_labels(true));
  ofs.close();

  CGAL::tetrahedral_isotropic_remeshing(t3,
    0.1,//Distance_from_corner_sizing_field(),
    CGAL::parameters::edge_is_constrained_map(Constraints_pmap(constraints))
                     .number_of_iterations(nb_iter)
                     .remesh_boundaries(true));

  std::ofstream ofs1("out1.mesh");
  CGAL::IO::write_MEDIT(ofs1, t3);
  ofs1.close();

  return EXIT_SUCCESS;
}
