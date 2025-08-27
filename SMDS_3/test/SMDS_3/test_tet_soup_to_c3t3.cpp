#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/tetrahedron_soup_to_triangulation_3.h>
#include <CGAL/IO/File_medit.h>

#include <vector>
#include <unordered_map>
#include <boost/unordered_map.hpp>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using DT3 = CGAL::Delaunay_triangulation_3<K>;
using Remeshing_triangulation = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;

using Point_3 = K::Point_3;
using Tetrahedron_3 = K::Tetrahedron_3;
using Vertex_handle = DT3::Vertex_handle;
using Subdomain_index = Remeshing_triangulation::Cell::Subdomain_index;
using Surface_patch_index = Remeshing_triangulation::Cell::Surface_patch_index;

int main(int , char* [])
{
  std::cout << "Random seed " << CGAL::get_default_random().get_seed() << std::endl;

  const int nbv = 100;

  //a triangulation
  DT3 delaunay;
  std::unordered_map<Vertex_handle, int> v2i;
  std::vector<DT3::Point>                points(nbv);
  std::vector<Tetrahedron_3>             tetrahedra;
  std::vector<std::array<int, 4> >       tets_by_indices;
  std::vector<Subdomain_index>           subdomains;
  boost::unordered_map<std::array<int, 3>, Surface_patch_index> border_facets;

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

  for (DT3::Cell_handle c : delaunay.all_cell_handles())
  {
    if(!delaunay.is_infinite(c))
      continue;

    //convex hull facet
    const int infinite_index = c->index(delaunay.infinite_vertex());
    DT3::Facet f(c, infinite_index);

    std::array<int, 3> facet;
    facet[0] = v2i.at(c->vertex((infinite_index + 1) % 4));
    facet[1] = v2i.at(c->vertex((infinite_index + 2) % 4));
    facet[2] = v2i.at(c->vertex((infinite_index + 3) % 4));
    border_facets[facet] = Surface_patch_index(12);
  }

  //build triangulation from tetrahedra
  Remeshing_triangulation tr
    = CGAL::tetrahedron_soup_to_triangulation_3<Remeshing_triangulation>(tetrahedra);
  assert(tr.is_valid());

  //build triangulation from indices
  Remeshing_triangulation tr2
    = CGAL::tetrahedron_soup_to_triangulation_3<Remeshing_triangulation>(
        points, tets_by_indices);
  assert(tr2.is_valid());

  //build triangulation from indices and subdomains
  Remeshing_triangulation tr3
    = CGAL::tetrahedron_soup_to_triangulation_3<Remeshing_triangulation>(
        points, tets_by_indices, CGAL::parameters::subdomain_indices(std::cref(subdomains)));
  assert(tr3.is_valid());

  //build triangulation from indices, subdomains and surface facets
  Remeshing_triangulation tr4
    = CGAL::tetrahedron_soup_to_triangulation_3<Remeshing_triangulation>(
        points, tets_by_indices,
        CGAL::parameters::subdomain_indices(std::cref(subdomains))
                         .surface_facets(std::cref(border_facets)));
  assert(tr4.is_valid());

  return EXIT_SUCCESS;
}
