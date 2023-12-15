#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/tetrahedron_soup_to_triangulation_3.h>
#include <CGAL/IO/File_medit.h>

#include <vector>
#include <unordered_map>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using DT3 = CGAL::Delaunay_triangulation_3<K>;
using Remeshing_triangulation = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;
using C3T3 = CGAL::Mesh_complex_3_in_triangulation_3<Remeshing_triangulation>;

using Point_3 = K::Point_3;
using Tetrahedron_3 = K::Tetrahedron_3;
using Vertex_handle = DT3::Vertex_handle;
using Subdomain_index = C3T3::Subdomain_index;

int main(int , char* [])
{
  const int nbv = 100;

  //a triangulation
  DT3 delaunay;
  std::unordered_map<Vertex_handle, int> v2i;
  std::vector<DT3::Point>                points(nbv);
  std::vector<Tetrahedron_3>             tetrahedra;
  std::vector<std::array<int, 4> >       tets_by_indices;
  std::vector< Subdomain_index>          subdomains;

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
  subdomains.reserve(delaunay.number_of_finite_cells());
  for (DT3::Cell_handle c : delaunay.finite_cell_handles())
  {
    tetrahedra.push_back(delaunay.tetrahedron(c));
    tets_by_indices.push_back( { v2i.at(c->vertex(0)),
                                 v2i.at(c->vertex(1)),
                                 v2i.at(c->vertex(2)),
                                 v2i.at(c->vertex(3)) } );
    subdomains.push_back(Subdomain_index(1));
  }

  //build triangulation from tetrahedra
  Remeshing_triangulation tr
    = CGAL::tetrahedron_soup_to_triangulation_3<Remeshing_triangulation>(tetrahedra);

  //build triangulation from indices
  Remeshing_triangulation tr2
    = CGAL::tetrahedron_soup_to_triangulation_3<Remeshing_triangulation>(
        points, tets_by_indices,
        CGAL::parameters::subdomain_indices(std::cref(subdomains)));

  //build a C3T3
  C3T3 c3t3;
  c3t3.triangulation() = tr;

  std::ofstream ofs("c3t3_output.mesh");
  CGAL::IO::write_MEDIT(ofs, c3t3);

  return EXIT_SUCCESS;
}
