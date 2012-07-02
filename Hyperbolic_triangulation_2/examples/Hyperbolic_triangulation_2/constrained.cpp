#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_hyperbolic_traits_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <cassert>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_hyperbolic_traits_2<K> Gt;

typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Gt> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
// typedef CGAL::No_intersection_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Gt, TDS> CDT;
typedef CDT::Point Point;

int main()
{
  CDT cdt;
  std::cout << "Inserting 8 constraints " << std::endl;
  
  // cdt.insert_constraint( Point(-0.25, 0), Point(0, 0.25) );
  // cdt.insert_constraint( Point( 0, 0.25), Point(0.25, 0) );
  // cdt.insert_constraint( Point( 0.25, 0), Point(0, -0.25) );
  // cdt.insert_constraint( Point( 0, -0.25), Point(-0.25, 0) );
  
  cdt.insert_constraint( Point(-0.85, 0), Point(-0.6, 0.6) );
  cdt.insert_constraint( Point( -0.6, 0.6), Point(0, 0.85) );
  cdt.insert_constraint( Point( 0, 0.85), Point(0.6, 0.6) );
  cdt.insert_constraint( Point( 0.6, 0.6), Point(0.85, 0) );
  cdt.insert_constraint( Point( 0.85, 0), Point(0.6, -0.6) );
  cdt.insert_constraint( Point( 0.6, -0.6), Point(0, -0.85) );
  cdt.insert_constraint( Point( 0, -0.85), Point(-0.6, -0.6) );
  cdt.insert_constraint( Point( -0.6, -0.6), Point(-0.85, 0) );
  
  //cdt.insert( Point( 0, 0 ) );
  
  // for (int i = 1; i < 6; ++i)
    // cdt.insert_constraint( Point(0,i), Point(6,i));
  // for (int j = 1; j < 6; ++j)
    // cdt.insert_constraint( Point(j,0), Point(j,6));

  assert(cdt.is_valid());
  int count = 0;
  int edgesCount = 0;
  for (CDT::Finite_edges_iterator eit = cdt.finite_edges_begin();
       eit != cdt.finite_edges_end();
       ++eit) {
    edgesCount++;
    if (cdt.is_constrained(*eit)) ++count;
  }
  std::cout << "The number of resulting constrained edges is  ";
  std::cout <<  count << std::endl;
  std::cout << "Edges count: ";
  std::cout << edgesCount << std::endl;
  
  int verticesCount = 0;
  for(CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
      vit != cdt.finite_vertices_end();
      ++vit) {
    std::cout << (*vit).point().x() << " " << (*vit).point().y() << std::endl;
    verticesCount++;
  }
  std::cout << verticesCount << std::endl;
  
  return 0;
}
