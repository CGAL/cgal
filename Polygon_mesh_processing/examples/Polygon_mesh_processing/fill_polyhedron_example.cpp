#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/triangulate_hole.h>

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/function_output_iterator.hpp>

struct Nop_functor 
{
  template<class T>
  void operator()(const T & /*t*/) const {}
};
typedef boost::function_output_iterator<Nop_functor> Nop_out;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>     Polyhedron;
typedef Polyhedron::Halfedge_iterator  Halfedge_iterator;
typedef Polyhedron::Facet_handle       Facet_handle;
typedef Polyhedron::Vertex_handle      Vertex_handle;

int main() {
  Polyhedron poly_1;
  std::ifstream input("data/max.off");
  if ( !input || !(input >> poly_1) || poly_1.empty() ) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }
  Polyhedron poly_2(poly_1), poly_3(poly_1);

  for(Halfedge_iterator h = poly_1.halfedges_begin(); h != poly_1.halfedges_end(); ++h) {
    if(h->is_border()) {
      std::vector<Facet_handle> patch;
      CGAL::Polygon_mesh_processing::triangulate_hole(poly_1,
        h, back_inserter(patch));
      std::cout << "Number of facets in constructed patch: " << patch.size() << std::endl;
    }
  }

  for(Halfedge_iterator h = poly_2.halfedges_begin(); h != poly_2.halfedges_end(); ++h) {
    if(h->is_border()) {
      std::vector<Facet_handle>  patch_facets;
      std::vector<Vertex_handle> patch_vertices;
      CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(poly_2,
                                        h,
                                        back_inserter(patch_facets),
                                        back_inserter(patch_vertices));
      std::cout << "Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << "Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
    }
  }

  for(Halfedge_iterator h = poly_3.halfedges_begin(); h != poly_3.halfedges_end(); ++h) {
    if(h->is_border()) {
      std::vector<Facet_handle>  patch_facets;
      std::vector<Vertex_handle> patch_vertices;
      bool success = CGAL::cpp11::get<0>(CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                                                            poly_3,
                                                            h,
                                                            back_inserter(patch_facets),
                                                            back_inserter(patch_vertices)));

      std::cout << "Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << "Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
      std::cout << "Is fairing successful: " << success << std::endl;
    }
  }
}
