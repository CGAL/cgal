#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/tags.h>
#include <utility>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Subdomain_index = int;
using Surface_patch_index = std::pair<int, int>;
using Curve_index = char;
using Corner_index = short;

using Cb = CGAL::Simplicial_mesh_cell_base_3<K, Subdomain_index, Surface_patch_index>;
using Vb = CGAL::Simplicial_mesh_vertex_base_3<K, Subdomain_index, Surface_patch_index,
                                                  Curve_index, Corner_index>;

using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag>;
using Tr = CGAL::Triangulation_3<K, Tds>;

using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

int main ()
{
  Tr triangulation;
  auto rng = CGAL::Random_points_on_sphere_3<K::Point_3>(1.);
  while (triangulation.number_of_vertices() < 100)
  {
    triangulation.insert(*rng++);
  }

  C3t3 c3t3;
  c3t3.triangulation() = std::move(triangulation);


  return EXIT_SUCCESS;
}
