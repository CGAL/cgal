#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#if CGAL_USE_CORE || CGAL_USE_LEDA
#  include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#endif
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>

#include <CGAL/IO/File_poly.h>

#include <fstream>
#include <iostream>
#include <cassert>

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

template <typename K>
struct Tester {
  void operator()() const {
    typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
    typedef typename CDT::Point Point;
    typedef typename CDT::size_type size_type;

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
      double x1,y1,x2,y2;
      edg_file >> x1 >> y1 >> x2 >> y2;
      Point p1(x1,y1), p2(x2,y2);
      cdt.insert_constraint(p1, p2);
    }

    const size_type number_of_vertices_edg = cdt.number_of_vertices();
    const size_type number_of_constrained_edges_edg = 
      number_of_constrained_edges(cdt);

    std::cout << "number of vertices: " << number_of_vertices_edg
              << "\nnumber of constrained edges: " 
              << number_of_constrained_edges_edg << "\n\n";


    // check that numbers of constrained edges and vertices are the same

    assert( number_of_constrained_edges_edg == 
                    number_of_constrained_edges_poly );

    assert( number_of_vertices_edg == number_of_vertices_poly );

    // CONFORMING

    std::cout << "Conforming...\n";
    CDT cdt2=cdt;
    assert(cdt.is_valid());

    // Gabriel for cdt
    make_conforming_Gabriel_2(cdt);
    assert(cdt.is_valid());
    CGAL::Triangulation_conformer_2<CDT> conformer(cdt);
    assert( conformer.is_conforming_Gabriel() );
    // Delaunay, then Gabriel for cdt2
    make_conforming_Delaunay_2(cdt2);
    CGAL::Triangulation_conformer_2<CDT> conformer2(cdt2);
    assert( conformer2.is_conforming_Delaunay() );    
    assert(cdt2.is_valid());
    std::cout << "Number of vertices after make_conforming_Delaunay_2: "
              << cdt2.number_of_vertices() << "\n";
    assert( cdt2.number_of_vertices() == 17 );
    assert(cdt2.is_valid());
    make_conforming_Gabriel_2(cdt2);
    std::cout << "Number of vertices after make_conforming_Gabriel_2: "
              << cdt2.number_of_vertices() << "\n";
    assert(cdt2.is_valid());
    assert( cdt2.number_of_vertices() == 45 );

    // check that numbers of vertices are the same in cdt and cdt2
    assert( cdt2.number_of_vertices() ==
                    cdt.number_of_vertices() );
  };
};


struct K_e_i : public CGAL::Exact_predicates_inexact_constructions_kernel {};
#if CGAL_USE_CORE || CGAL_USE_LEDA
struct K_e_e : public CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt {};
#endif

int main()
{
  std::cerr << "TESTING WITH Exact_predicates_inexact_constructions_kernel...\n\n";
  Tester<K_e_i> tester;
  tester();
#if CGAL_USE_CORE || CGAL_USE_LEDA
  std::cerr << "\n\nTESTING WITH Exact_predicates_exact_constructions_kernel_with_sqrt...\n\n";
  Tester<K_e_e> tester2;
  tester2();
#endif
}
