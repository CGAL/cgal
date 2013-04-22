#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_triangulation_traits_2<K> GT;

typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT> PDT;

typedef PDT::Face_handle       Face_handle;
typedef PDT::Vertex_handle     Vertex_handle;
typedef PDT::Locate_type       Locate_type;
typedef PDT::Point             Point;
typedef PDT::Iso_rectangle     Iso_rectangle;

int main()
{
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  Iso_rectangle domain(-1, -1, 2, 2); // The cube for the periodic domain

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  // construction from a list of points :
  std::list<Point> L;
  L.push_front(Point(0, 0));
  L.push_front(Point(1, 0));
  L.push_front(Point(0, 1));

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  PDT T(L.begin(), L.end(), domain); // Put the domain with the constructor

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  int n = T.number_of_vertices();

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  // insertion from a vector :
  std::vector<Point> V(3);
  V[0] = Point(0, 0);
  V[1] = Point(1, 1);
  V[2] = Point(-1, -1);

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  n = n + T.insert(V.begin(), V.end());

  assert( n == 5 );       // 6 points have been inserted, one is a duplicate
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  assert( T.is_valid() ); // checking validity of T

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  Locate_type lt;
  int li;
  Point p(0, 0);
  Face_handle fh = T.locate(p, lt, li);
  // p is the vertex of c of index li :
  assert( lt == PDT::VERTEX );
  assert( fh->vertex(li)->point() == p );

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  Vertex_handle v = fh->vertex( (li + 1) & 3 );
  // v is another vertex of c
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  Face_handle nb = fh->neighbor(li);
  // nb = neighbor of fh opposite to the vertex associated with p
  // nb must have vertex v :
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  int nli;
  assert( nb->has_vertex( v, nli ) );
  // nli is the index of v in nc

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  std::ofstream oFileT("output.tri", std::ios::out);
  // writing file output;
  oFileT << T;

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  PDT T1;
  std::ifstream iFileT("output.tri", std::ios::in);
  // reading file output;
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  iFileT >> T1;
  assert( T1.is_valid() );
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  assert( T1.number_of_vertices() == T.number_of_vertices() );
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  assert( T1.number_of_faces() == T.number_of_faces() );
  std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  return 0;
}
