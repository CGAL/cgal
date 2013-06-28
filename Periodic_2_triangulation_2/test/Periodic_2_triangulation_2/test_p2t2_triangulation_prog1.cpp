// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_triangulation_2.h>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Periodic_2_triangulation_traits_2<K>  Gt;

typedef CGAL::Periodic_2_triangulation_2<Gt>        Triangulation;
typedef Triangulation::Vertex_circulator            Vertex_circulator;
typedef Triangulation::Point                        Point;
typedef Gt::Vector_2                                Vector;
typedef Triangulation::Segment                      Segment;
typedef Triangulation::Triangle                     Triangle;
typedef Triangulation::Periodic_point_iterator      Periodic_point_iterator;
typedef Triangulation::Periodic_segment_iterator    Periodic_segment_iterator;
typedef Triangulation::Periodic_triangle_iterator   Periodic_triangle_iterator;

void test_point_location_in_a_triangulation_with_a_single_point(Triangulation &t)
{
  CGAL_assertion(t.number_of_vertices() == 1);

  Point &p = t.vertices_begin()->point();
  Triangulation::Locate_type lt;
  int li;

  for (int offset_x = 0; offset_x < t.number_of_sheets()[0]; ++offset_x)
    {
      for (int offset_y = 0; offset_y < t.number_of_sheets()[1]; ++offset_y)
        {
          Triangulation::Offset offset(offset_x, offset_y);
          {
            // Test the triangle iterator

            Triangulation::Periodic_triangle_iterator triang_it;
            const Triangulation::Periodic_triangle_iterator triang_it_beyond = t.periodic_triangles_end();
            {
              // Locate the point at the only vertex
              int on_boundary = 0;
              int on_inside = 0;
              triang_it = t.periodic_triangles_begin();
              while(triang_it != triang_it_beyond)
                {
                  CGAL::Bounded_side side = t.side_of_face(p, offset, triang_it.get_face(), lt, li);
                  if (side != CGAL::ON_UNBOUNDED_SIDE)
                    {
                      CGAL_assertion(lt == Triangulation::VERTEX);
                    }
                  on_boundary += (side == CGAL::ON_BOUNDARY ? 1 : 0);
                  on_inside += (side == CGAL::ON_BOUNDED_SIDE ? 1 : 0);
                  Triangle triangle = t.triangle(*triang_it);
                  triangle.vertex(0); // Avoid warning
                  ++triang_it;
                }
              CGAL_assertion(on_boundary == 6);
              CGAL_assertion(on_inside == 0);
            }
            {
              // Locate point on an edge
              int on_boundary = 0;
              int on_inside = 0;
              triang_it = t.periodic_triangles_begin();
              while(triang_it != triang_it_beyond)
                {
                  CGAL::Bounded_side side = t.side_of_face(p + Vector(0.0, 0.1), offset, triang_it.get_face(), lt, li);
                  if (side != CGAL::ON_UNBOUNDED_SIDE)
                    {
                      CGAL_assertion(lt == Triangulation::EDGE);
                    }
                  on_boundary += (side == CGAL::ON_BOUNDARY ? 1 : 0);
                  on_inside += (side == CGAL::ON_BOUNDED_SIDE ? 1 : 0);
                  Triangle triangle = t.triangle(*triang_it);
                  triangle.vertex(0); // Avoid warning
                  ++triang_it;
                }
              CGAL_assertion(on_boundary == 2);
              CGAL_assertion(on_inside == 0);
            }
            {
              // Locate point inside a face
              int on_boundary = 0;
              int on_inside = 0;
              triang_it = t.periodic_triangles_begin();
              while(triang_it != triang_it_beyond)
                {
                  CGAL::Bounded_side side = t.side_of_face(p + Vector(0.1, 0.2), offset, triang_it.get_face(), lt, li);
                  if (side != CGAL::ON_UNBOUNDED_SIDE)
                    {
                      CGAL_assertion(lt == Triangulation::FACE);
                    }
                  on_boundary += (side == CGAL::ON_BOUNDARY ? 1 : 0);
                  on_inside += (side == CGAL::ON_BOUNDED_SIDE ? 1 : 0);
                  Triangle triangle = t.triangle(*triang_it);
                  triangle.vertex(0); // Avoid warning
                  ++triang_it;
                }
              CGAL_assertion(on_boundary == 0);
              CGAL_assertion(on_inside == 1);
            }
          }
          {
            // Test the face iterator
            Triangulation::Face_iterator face_it;
            const Triangulation::Face_iterator face_it_beyond = t.faces_end();

            {
              // Locate the point at the only vertex
              int on_boundary = 0;
              int on_inside = 0;
              face_it = t.faces_begin();
              while(face_it != face_it_beyond)
                {
                  CGAL::Bounded_side side = t.side_of_face(p, offset, face_it, lt, li);
                  if (side != CGAL::ON_UNBOUNDED_SIDE)
                    {
                      CGAL_assertion(lt == Triangulation::VERTEX);
                    }
                  on_boundary += (side == CGAL::ON_BOUNDARY ? 1 : 0);
                  on_inside += (side == CGAL::ON_BOUNDED_SIDE ? 1 : 0);
                  ++face_it;
                }
              CGAL_assertion(on_boundary == 6);
              CGAL_assertion(on_inside == 0);
            }
            {
              // Locate point on an edge
              int on_boundary = 0;
              int on_inside = 0;
              face_it = t.faces_begin();
              while(face_it != face_it_beyond)
                {
                  CGAL::Bounded_side side = t.side_of_face(p + Vector(0.1, 0.0), offset, face_it, lt, li);
                  if (side != CGAL::ON_UNBOUNDED_SIDE)
                    {
                      CGAL_assertion(lt == Triangulation::EDGE);
                    }
                  on_boundary += (side == CGAL::ON_BOUNDARY ? 1 : 0);
                  on_inside += (side == CGAL::ON_BOUNDED_SIDE ? 1 : 0);
                  ++face_it;
                }
              CGAL_assertion(on_boundary == 2);
              CGAL_assertion(on_inside == 0);
            }
            {
              // Locate point inside a face
              int on_boundary = 0;
              int on_inside = 0;
              face_it = t.faces_begin();
              while(face_it != face_it_beyond)
                {
                  CGAL::Bounded_side side = t.side_of_face(p + Vector(0.1, 0.2), offset, face_it, lt, li);
                  if (side != CGAL::ON_UNBOUNDED_SIDE)
                    {
                      CGAL_assertion(lt == Triangulation::FACE);
                    }
                  on_boundary += (side == CGAL::ON_BOUNDARY ? 1 : 0);
                  on_inside += (side == CGAL::ON_BOUNDED_SIDE ? 1 : 0);
                  ++face_it;
                }
              CGAL_assertion(on_boundary == 0);
              CGAL_assertion(on_inside == 1);
            }
          }
        }
    }
}

int main()
{
  Triangulation t;

  // Insert the first point
  Point first_point(0.5, 0.5);
  t.insert(first_point);
  CGAL_assertion(t.is_valid(true));

  {
    // Testing the point iterator
    Periodic_point_iterator it = t.periodic_points_begin();
    Periodic_point_iterator beyond = t.periodic_points_end();
    CGAL_assertion(std::distance(it, beyond) == 9);
    while(it != beyond)
      {
        Point p = t.point(*it);
        p.x(); // Avoid warning
        ++it;
      }
  }
  {
    // Testing the segment iterator
    Periodic_segment_iterator it = t.periodic_segments_begin();
    Periodic_segment_iterator beyond = t.periodic_segments_end();
    CGAL_assertion(std::distance(it, beyond) == 27);
    while(it != beyond)
      {
        Segment s = t.segment(*it);
        s.point(0); // Avoid warning
        ++it;
      }
  }
  {
    // Testing the triangle iterator
    Periodic_triangle_iterator it = t.periodic_triangles_begin();
    Periodic_triangle_iterator beyond = t.periodic_triangles_end();
    CGAL_assertion(std::distance(it, beyond) == 18);
    while(it != beyond)
      {
        Triangle triangle = t.triangle(*it);
        triangle.vertex(0); // Avoid warning
        ++it;
      }
  }

  test_point_location_in_a_triangulation_with_a_single_point(t);

  //  std::ifstream in("data/triangulation_prog1.cin");
  //  std::istream_iterator<Point> begin(in);
  //  std::istream_iterator<Point> end;
  //  t.insert(begin, end);

  return 0;
}
