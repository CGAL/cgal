#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>

#include <CGAL/IO/File_poly.h>

#include <fstream>
#include <iostream>

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Point Point;
typedef CDT::Vertex_handle Vertex_handle;

typedef CDT::size_type size_type;

template <class CTr>
typename CTr::size_type number_of_constrained_edges(const CTr& tr)
{
  typename CTr::size_type nedges = 0;
  for(typename CTr::Finite_edges_iterator eit = tr.finite_edges_begin();
      eit != tr.finite_edges_end();
      ++eit)
    if(tr.is_constrained(*eit))
      ++nedges;
  return nedges;
}

int main()
{
  CDT cdt;

  // CHECK FIRST read_triangle_poly_file AND operator>> OF CDT

  // read a poly file

  std::cout << "Reading fish.poly...\n";
  std::ifstream poly_file("fish.poly");
  CGAL::read_triangle_poly_file(cdt, poly_file);

  const size_type number_of_vertices_poly = cdt.number_of_vertices();
  const size_type number_of_constrained_edges_poly = 
    number_of_constrained_edges(cdt);

  std::cout << "number of vertices: " << number_of_vertices_poly
            << "\nnumber of constrained edges: " 
            << number_of_constrained_edges_poly << "\n\n";

  // read a CGAL file (edg file).

  std::cout << "Reading fish.edg...\n";

  std::ifstream edg_file("fish.edg");

  cdt.clear();

  size_type nedges = 0;
  edg_file >> nedges;
  for(size_type n = 0; n < nedges; ++n) {
    Point p1, p2;
    edg_file >> p1 >> p2;
    cdt.insert_constraint(p1, p2);
  }

  const size_type number_of_vertices_edg = cdt.number_of_vertices();
  const size_type number_of_constrained_edges_edg = 
    number_of_constrained_edges(cdt);

  std::cout << "number of vertices: " << number_of_vertices_edg
            << "\nnumber of constrained edges: " 
            << number_of_constrained_edges_edg << "\n\n";


  // check that numbers of constrained edges and vertices are the same

  CGAL_assertion( number_of_constrained_edges_edg == 
                  number_of_constrained_edges_poly );

  CGAL_assertion( number_of_vertices_edg == number_of_vertices_poly );

  // CONFORMING

  std::cout << "Conforming...\n";
  CDT cdt2=cdt;

  // Gabriel for cdt
  make_conforming_Gabriel_2(cdt);
  // Delaunay, then Gabriel for cdt2
  make_conforming_Delaunay_2(cdt2);
  std::cout << "Number of vertices after make_conforming_Delaunay_2: "
            << cdt2.number_of_vertices() << "\n";
  CGAL_assertion( cdt2.number_of_vertices() == 17 );
  make_conforming_Gabriel_2(cdt2);
  std::cout << "Number of vertices after make_conforming_Gabriel_2: "
            << cdt2.number_of_vertices() << "\n";
  CGAL_assertion( cdt2.number_of_vertices() == 45 );

  // check that numbers of vertices are the same in cdt and cdt2
  CGAL_assertion( cdt2.number_of_vertices() ==
                  cdt.number_of_vertices() );

}
