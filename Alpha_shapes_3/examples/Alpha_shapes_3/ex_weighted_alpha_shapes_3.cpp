#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Regular_triangulation_3.h>

#include <iostream>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_vertex_base_3<K>        Rvb;
typedef CGAL::Alpha_shape_vertex_base_3<K,Rvb>              Vb;

typedef CGAL::Regular_triangulation_cell_base_3<K>          Rcb;
typedef CGAL::Alpha_shape_cell_base_3<K,Rcb>                Cb;

typedef CGAL::Triangulation_data_structure_3<Vb,Cb>         Tds;
typedef CGAL::Regular_triangulation_3<K,Tds>                Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>                Alpha_shape_3;

typedef Alpha_shape_3::Cell_handle                          Cell_handle;
typedef Alpha_shape_3::Vertex_handle                        Vertex_handle;
typedef Alpha_shape_3::Facet                                Facet;
typedef Alpha_shape_3::Edge                                 Edge;
typedef Triangulation_3::Weighted_point                     Weighted_point;
typedef Triangulation_3::Bare_point                         Bare_point;

int main()
{
  std::list<Weighted_point> lwp;

  //input : a small molecule
  lwp.push_back(Weighted_point(Bare_point( 1, -1, -1), 4));
  lwp.push_back(Weighted_point(Bare_point(-1,  1, -1), 4));
  lwp.push_back(Weighted_point(Bare_point(-1, -1,  1), 4));
  lwp.push_back(Weighted_point(Bare_point( 1,  1,  1), 4));
  lwp.push_back(Weighted_point(Bare_point( 2,  2,  2), 1));

  //build alpha_shape  in GENERAL mode and set alpha=0
  Alpha_shape_3  as(lwp.begin(), lwp.end(), 0, Alpha_shape_3::GENERAL);

  //explore the 0-shape - It is dual to the boundary of the union.
  std::list<Cell_handle> cells;
  std::list<Facet>       facets;
  std::list<Edge>        edges;
  as.get_alpha_shape_cells(std::back_inserter(cells),
                           Alpha_shape_3::INTERIOR);
  as.get_alpha_shape_facets(std::back_inserter(facets),
                            Alpha_shape_3::REGULAR);
  as.get_alpha_shape_facets(std::back_inserter(facets),
                            Alpha_shape_3::SINGULAR);
  as.get_alpha_shape_edges(std::back_inserter(edges),
                           Alpha_shape_3::SINGULAR);
  std::cout << " The 0-shape has : " << std::endl;
  std::cout << cells.size() << " interior tetrahedra" << std::endl;
  std::cout << facets.size() << " boundary facets" << std::endl;
  std::cout << edges.size()  << " singular edges" << std::endl;
  return 0;
}
