// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include "./types.h"

Face_handle test_point_location(const Triangulation &t,
                                const Point &query,
                                const Triangulation::Locate_type &lt_in)
{
  Triangulation::Locate_type lt, lt2;
  int li, li2;
  Face_handle fh;
  CGAL::Bounded_side bs;
  CGAL::Oriented_side os;

  fh = t.locate(query, lt, li);
  assert(lt == lt_in);
  if (lt_in == Triangulation::EMPTY) {
    assert(fh == Face_handle());
    return fh;
  }

  bs = t.side_of_face(query, fh, lt2, li2);
  os = t.oriented_side(fh, query);
  CGAL_USE(bs);
  CGAL_USE(os);

  assert(lt2 == lt_in);

  switch (lt_in)
  {
    case Triangulation::VERTEX:
    case Triangulation::EDGE:
    {
      assert(fh != Face_handle());
      assert(bs == CGAL::ON_BOUNDARY);
      assert(os == CGAL::ON_ORIENTED_BOUNDARY);

      assert(li == li2);
      break;
    }
    case Triangulation::FACE:
    {
      assert(fh != Face_handle());
      assert(bs == CGAL::ON_BOUNDED_SIDE);
      assert(os == CGAL::ON_POSITIVE_SIDE);
      break;
    }
    default:
    {
      // Handled above
      assert(false);
      break;
    }
  }

  return fh;
}

int main()
{
  Triangulation t;
  Face_handle fh;

  // Check the empty triangulation
  fh = test_point_location(t, Point(0.5, 0.5), Triangulation::EMPTY);
  assert(fh == Face_handle());

  // Insert the first point
  Point p0(0.5, 0.5);
  Vertex_handle vh0 = t.insert(p0);
  assert(t.is_valid(true));
  CGAL_USE(vh0);

  fh = test_point_location(t, p0, Triangulation::VERTEX);
  assert(fh->has_vertex(vh0));

  fh = test_point_location(t, p0 + Vector(0.1, 0.1), Triangulation::EDGE);
  assert(fh->has_vertex(vh0));

  fh = test_point_location(t, p0 + Vector(-0.1, -0.1), Triangulation::EDGE);
  assert(fh->has_vertex(vh0));

  fh = test_point_location(t, p0 + Vector(-0.2, -0.3), Triangulation::FACE);
  assert(fh->has_vertex(vh0));

  assert(t.is_valid(true));

  // Insert the second point on an edge
  Point p1(0.7, 0.7);
  Vertex_handle vh1 = t.insert(p1);
  CGAL_USE(vh1);
  assert(t.is_valid(true));

  fh = test_point_location(t, p0, Triangulation::VERTEX);
  assert(fh->has_vertex(vh0));

  fh = test_point_location(t, p1, Triangulation::VERTEX);
  assert(fh->has_vertex(vh1));

  fh = test_point_location(t, p0 + Vector(0.1, 0.1), Triangulation::EDGE);
  assert(fh->has_vertex(vh0));
  assert(fh->has_vertex(vh1));

  fh = test_point_location(t, p0 + Vector(-0.1, -0.1), Triangulation::EDGE);
  assert(fh->has_vertex(vh0));
  assert(!fh->has_vertex(vh1));

  fh = test_point_location(t, p1 + Vector(0.1, 0.1), Triangulation::EDGE);
  assert(!fh->has_vertex(vh0));
  assert(fh->has_vertex(vh1));

  fh = test_point_location(t, p0 + Vector(-0.02, -0.03), Triangulation::FACE);
  assert(fh->has_vertex(vh0));

  fh = test_point_location(t, p1 + Vector(-0.02, -0.03), Triangulation::FACE);
  assert(fh->has_vertex(vh1));

  assert(t.is_valid(true));

  // Insert the third point in a face
  Point p2(0.8, 0.6);
  Vertex_handle vh2 = t.insert(p2);
  CGAL_USE(vh2);
  assert(t.is_valid(true));

  fh = test_point_location(t, p0, Triangulation::VERTEX);
  assert(fh->has_vertex(vh0));
  fh = test_point_location(t, p1, Triangulation::VERTEX);
  assert(fh->has_vertex(vh1));
  fh = test_point_location(t, p2, Triangulation::VERTEX);
  assert(fh->has_vertex(vh2));

  fh = test_point_location(t, Point(0.6, 0.6), Triangulation::EDGE);
  assert(fh->has_vertex(vh0));
  assert(fh->has_vertex(vh1));

  test_point_location(t, Point(0.7, 0.6), Triangulation::FACE);
  test_point_location(t, p0 + Vector(-0.02, -0.03), Triangulation::FACE);
  test_point_location(t, p0 + Vector(0.02, -0.03), Triangulation::FACE);
  test_point_location(t, p0 + Vector(-0.02, 0.03), Triangulation::FACE);
  test_point_location(t, p0 + Vector(0.02, 0.03), Triangulation::FACE);

  return 0;
}
