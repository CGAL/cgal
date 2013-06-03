// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include "./types.h"

int main()
{
  Point p;
  Triangulation t;

  CGAL_assertion(t.is_valid());

  Point p0(0.5, 0.5);
  CGAL_assertion_code(Vertex_handle vh0 = )
  t.insert(p0);
  CGAL_assertion(t.is_valid());

  CGAL_assertion(t.number_of_vertices() == 1);
  CGAL_assertion(vh0 == t.insert(p0));
  CGAL_assertion(t.is_valid());

  Point p1(0.6, 0.5);
  CGAL_assertion_code(Vertex_handle vh1 = )
  t.insert(p1);
  CGAL_assertion(t.is_valid());
  CGAL_assertion(t.number_of_vertices() == 2);
  CGAL_assertion(vh0 != vh1);
  CGAL_assertion(vh0 == t.insert(p0));
  CGAL_assertion(vh1 == t.insert(p1));
  CGAL_assertion(t.is_valid());

  Point p2(0.5, 0.7);
  CGAL_assertion_code(Vertex_handle vh2 = )
  t.insert(p2);
  CGAL_assertion(t.is_valid());
  CGAL_assertion(t.number_of_vertices() == 3);
  CGAL_assertion(vh0 == t.insert(p0));
  CGAL_assertion(vh1 == t.insert(p1));
  CGAL_assertion(vh2 == t.insert(p2));
  CGAL_assertion(t.is_valid());

  Point p3(0.3, 0.4);
  CGAL_assertion_code(Vertex_handle vh3 = )
  t.insert(p3);
  CGAL_assertion(t.is_valid());
  CGAL_assertion(t.number_of_vertices() == 4);
  CGAL_assertion(vh0 == t.insert(p0));
  CGAL_assertion(vh1 == t.insert(p1));
  CGAL_assertion(vh2 == t.insert(p2));
  CGAL_assertion(vh3 == t.insert(p3));
  CGAL_assertion(t.is_valid());

  return 0;
}
