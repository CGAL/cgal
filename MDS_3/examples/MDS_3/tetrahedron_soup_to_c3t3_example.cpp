#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/tetrahedron_soup_to_triangulation_3.h>

#include <vector>
#include <boost/unordered_map.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_3<K> DT3;
typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Remeshing_triangulation> C3T3;

typedef K::Point_3         Point_3;
typedef K::Tetrahedron_3   Tetrahedron_3;
typedef DT3::Vertex_handle Vertex_handle;

int main(int , char* [])
{
  const int nbv = 100;

  //a triangulation
  DT3 delaunay;
  boost::unordered_map<Vertex_handle, int>     v2i;
  std::vector<DT3::Point>                points(nbv);
  std::vector<Tetrahedron_3>             tetrahedra;
  std::vector<std::array<int, 5> >       tets_by_indices;

  //insert random points
  CGAL::Random_points_in_cube_3<Point_3> randp(2.);
  int i = 0;
  while (i < nbv)
  {
    points[i] = *randp++;
    Vertex_handle v = delaunay.insert(points[i]);
    v2i[v] = i++;
  }

  tetrahedra.reserve(delaunay.number_of_finite_cells());
  tets_by_indices.reserve(delaunay.number_of_finite_cells());
  for (DT3::Cell_handle c : delaunay.finite_cell_handles())
  {
    tetrahedra.push_back(delaunay.tetrahedron(c));

    std::array<int, 5> tet;
    tet[0] = v2i.at(c->vertex(0));
    tet[1] = v2i.at(c->vertex(1));
    tet[2] = v2i.at(c->vertex(2));
    tet[3] = v2i.at(c->vertex(3));
    tet[4] = Remeshing_triangulation::Cell::Subdomain_index(1);

    tets_by_indices.push_back(tet);
  }

  //build triangulation from tetrahedra
  Remeshing_triangulation tr;
  CGAL::tetrahedron_soup_to_triangulation_3(tetrahedra, tr);

  //buid triangulation from indices
  Remeshing_triangulation tr2;
  CGAL::tetrahedron_soup_to_triangulation_3(points, tets_by_indices, tr2);

  //build a C3T3
  C3T3 c3t3;
  c3t3.triangulation() = tr;

  std::ofstream ofs("c3t3_output.mesh");
  c3t3.output_to_medit(ofs);

  return EXIT_SUCCESS;
}
