// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include "./types.h"

void insert_in_edge(Triangulation &t, const Point &p)
{
  Triangulation::Locate_type lt;
  int li;

  Face_handle fh = t.locate(p, lt, li);
  CGAL_assertion(lt == Triangulation::EDGE);
  t.insert(p, fh);
  CGAL_assertion(t.is_valid());
}

int main()
{
  Point p;
  Triangulation t;

  // Insert the first point
  t.insert(Point(0.5, 0.5));

  insert_in_edge(t, Point(0.5, 0.7));
  insert_in_edge(t, Point(0.7, 0.5));
  insert_in_edge(t, Point(0.7, 0.7));
  insert_in_edge(t, Point(0.8, 0.8));
  insert_in_edge(t, Point(0.6, 0.6));

  insert_in_edge(t, Point(0.3, 0.5));
  insert_in_edge(t, Point(0.5, 0.3));
  insert_in_edge(t, Point(0.3, 0.3));
  insert_in_edge(t, Point(0.4, 0.4));
  insert_in_edge(t, Point(0.2, 0.2));

  return 0;
}
