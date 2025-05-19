// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include "./types.h"

int main()
{
  Point p;
  Triangulation t;

  assert(t.is_valid());

  Point p0(0.5, 0.5);
  Vertex_handle vh0 = t.insert(p0);
  assert(t.is_valid());

  assert(t.number_of_vertices() == 1);
  assert(vh0 == t.insert(p0));
  assert(t.is_valid());

  Point p1(0.6, 0.5);
  Vertex_handle vh1 = t.insert(p1);
  assert(t.is_valid());
  assert(t.number_of_vertices() == 2);
  assert(vh0 != vh1);
  assert(vh0 == t.insert(p0));
  assert(vh1 == t.insert(p1));
  assert(t.is_valid());

  Point p2(0.5, 0.7);
  Vertex_handle vh2 = t.insert(p2);
  assert(t.is_valid());
  assert(t.number_of_vertices() == 3);
  assert(vh0 == t.insert(p0));
  assert(vh1 == t.insert(p1));
  assert(vh2 == t.insert(p2));
  assert(t.is_valid());

  Point p3(0.3, 0.4);
  Vertex_handle vh3 = t.insert(p3);
  assert(t.is_valid());
  assert(t.number_of_vertices() == 4);
  assert(vh0 == t.insert(p0));
  assert(vh1 == t.insert(p1));
  assert(vh2 == t.insert(p2));
  assert(vh3 == t.insert(p3));
  assert(t.is_valid());

  return 0;
}
