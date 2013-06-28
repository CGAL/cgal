#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>

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
  Iso_rectangle domain(-1, -1, 2, 2); // The cube for the periodic domain

  // construction from a list of points :
  std::list<Point> L;
  L.push_front(Point(0, 0));
  L.push_front(Point(1, 0));
  L.push_front(Point(0, 1));

  PDT T(L.begin(), L.end(), domain); // Put the domain with the constructor

  size_t n = T.number_of_vertices();

  // insertion from a vector :
  std::vector<Point> V(3);
  V[0] = Point(0, 0);
  V[1] = Point(1, 1);
  V[2] = Point(-1, -1);

  n = n + T.insert(V.begin(), V.end());

  assert( n == 5 );       // 6 points have been inserted, one is a duplicate
  assert( T.is_valid() ); // checking validity of T

  Locate_type lt;
  int li;
  Point p(0, 0);
  Face_handle fh = T.locate(p, lt, li);
  // p is the vertex of c of index li :
  assert( lt == PDT::VERTEX );
  assert( fh->vertex(li)->point() == p );

  Vertex_handle v = fh->vertex( (li + 1) % 3 );
  // v is another vertex of c
  Face_handle nb = fh->neighbor(li);
  // nb = neighbor of fh opposite to the vertex associated with p
  // nb must have vertex v :
  int nli;
  assert( nb->has_vertex( v, nli ) );
  // nli is the index of v in nc

  std::ofstream oFileT("output.tri", std::ios::out);
  // writing file output;
  oFileT << T;

  PDT T1;
  std::ifstream iFileT("output.tri", std::ios::in);
  // reading file output;
  iFileT >> T1;
  assert( T1.is_valid() );
  assert( T1.number_of_vertices() == T.number_of_vertices() );
  assert( T1.number_of_faces() == T.number_of_faces() );

  return 0;
}
