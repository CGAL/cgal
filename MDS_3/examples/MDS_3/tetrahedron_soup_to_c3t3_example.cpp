#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/tetrahedron_soup_to_triangulation_3.h>

#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_3<K> DT3;
typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

typedef K::Point_3       Point_3;
typedef K::Tetrahedron_3 Tetrahedron_3;


int main(int argc, char* argv[])
{
  //build a triangulation
  DT3 delaunay;
  CGAL::Random_points_in_cube_3<Point_3> randp(2.);
  while (delaunay.number_of_vertices() < 100)
    delaunay.insert(*randp++);

  //collect tetrahedra
  std::vector<Tetrahedron_3> tetrahedra(delaunay.number_of_finite_cells());
  for (DT3::Cell_handle c : delaunay.finite_cell_handles())
    tetrahedra.push_back(delaunay.tetrahedron(c));

  //build triangulation
  Remeshing_triangulation tr;
  CGAL::tetrahedron_soup_to_triangulation_3(tetrahedra, tr);

  std::ofstream os("dt_rebuilt.mesh");
  CGAL::write_MEDIT(os, tr);

  return EXIT_SUCCESS;
}
