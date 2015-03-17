#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/triangulate_hole.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/foreach.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>     Polyhedron;

typedef Polyhedron::Halfedge_handle    Halfedge_handle;
typedef Polyhedron::Facet_handle       Facet_handle;
typedef Polyhedron::Vertex_handle      Vertex_handle;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/mech-holes-shark.off";
  std::ifstream input(filename);

  Polyhedron poly_1;
  if ( !input || !(input >> poly_1) || poly_1.empty() ) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }
  Polyhedron poly_2(poly_1), poly_3(poly_1);

  // in poly_1, incrementally fill the holes
  BOOST_FOREACH(Halfedge_handle h, halfedges(poly_1))
  {
    if(h->is_border())
    {
      std::vector<Facet_handle> patch;
      CGAL::Polygon_mesh_processing::triangulate_hole(poly_1,
                                        h,
                                        std::back_inserter(patch));
      std::cout << "Number of facets in constructed patch: " << patch.size() << std::endl;
    }
  }

  // in poly_2, incrementally fill and refine the holes
  BOOST_FOREACH(Halfedge_handle h, halfedges(poly_2))
  {
    if(h->is_border())
    {
      std::vector<Facet_handle>  patch_facets;
      std::vector<Vertex_handle> patch_vertices;
      CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(poly_2,
                                        h,
                                        std::back_inserter(patch_facets),
                                        std::back_inserter(patch_vertices));
      std::cout << "Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << "Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
    }
  }

  // in poly_3, incrementally fill, refine and fair the holes
  BOOST_FOREACH(Halfedge_handle h, halfedges(poly_3))
  {
    if (h->is_border())
    {
      std::vector<Facet_handle>  patch_facets;
      std::vector<Vertex_handle> patch_vertices;
      bool success = CGAL::cpp11::get<0>(CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                                        poly_3,
                                        h,
                                        std::back_inserter(patch_facets),
                                        std::back_inserter(patch_vertices)));

      std::cout << "Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << "Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
      std::cout << "Is fairing successful: " << success << std::endl;
    }
  }

  return 0;
}
