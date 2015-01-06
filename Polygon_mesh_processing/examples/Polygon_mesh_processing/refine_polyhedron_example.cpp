#include <CGAL/Hole_filling.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/iterator.h>

#include <iostream>
#include <fstream>
#include <functional>

#include <boost/iterator/transform_iterator.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
typedef Polyhedron::Facet_handle    Facet_handle;
typedef Polyhedron::Vertex_handle   Vertex_handle;
typedef Polyhedron::Facet           Facet;

struct Facet_to_facet_handle 
  : public std::unary_function<Facet&, Facet_handle>
{
  result_type operator()(argument_type f) const
  { return f.halfedge()->facet(); }
};


int main() {
  Polyhedron poly_1;
  std::ifstream input("data/max.off");
  if ( !input || !(input >> poly_1) || poly_1.empty() ) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }
  Polyhedron poly_2 = poly_1;

  std::vector<Facet_handle>  new_facets;
  std::vector<Vertex_handle> new_vertices;
  // `facets_begin()` returns `Facet_iterator` which is an iterator over `Facet`, 
  // what `refine()` requires is an iterator over `Facet_handle`, hence a transformer is used
  CGAL::refine(poly_1, 
    boost::make_transform_iterator(poly_1.facets_begin(), Facet_to_facet_handle()),
    boost::make_transform_iterator(poly_1.facets_end()  , Facet_to_facet_handle()),
    back_inserter(new_facets), back_inserter(new_vertices), 2.0);
  
  std::ofstream poly_1_off("data/poly_1.off");
  poly_1_off << poly_1;
  poly_1_off.close();
  
  CGAL::refine(poly_2, 
    boost::make_transform_iterator(poly_2.facets_begin(), Facet_to_facet_handle()),
    boost::make_transform_iterator(poly_2.facets_end()  , Facet_to_facet_handle()),
               CGAL::Emptyset_iterator(), CGAL::Emptyset_iterator(), 3.0);
	
  std::ofstream poly_2_off("data/poly_2.off");
  poly_2_off << poly_2;
  poly_2_off.close();
}
