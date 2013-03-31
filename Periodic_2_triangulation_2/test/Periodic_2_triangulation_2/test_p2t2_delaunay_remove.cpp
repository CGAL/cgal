// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#define CGAL_USE_DELAUNAY 1

#include "./types.h"

const int N_PTS = 75;

int main()
{
  Triangulation t;

  Vector midpoint(0.5, 0.5);

  Face_handle fh;
  Triangulation::Locate_type lt;
  int i;
  fh = t.locate(Point(0, 0) + midpoint, lt, i);
  CGAL_assertion(lt == Triangulation::EMPTY);

  Vertex_handle vh_midpoint = t.insert(Point(0, 0) + midpoint);
  fh = t.locate(Point(0, 0) + midpoint, lt, i);
  CGAL_assertion(lt == Triangulation::VERTEX && fh->vertex(i) == vh_midpoint);
  t.remove(vh_midpoint);
  CGAL_assertion(t.empty());

  // High degree vertex
  for (int n = 3; n < 8; ++n)
    {
      vh_midpoint = t.insert(Point(0, 0) + midpoint);
      for (int i = 0; i < n; ++i)
        {
          t.insert(Point(0.3 * sin(i * 1.0 / n * 2 * M_PI), 0.3 * cos(i * 1.0 / n * 2 * M_PI)) + midpoint);
        }
      t.remove(vh_midpoint);
      CGAL_assertion(t.is_valid(true));
      while (!t.empty())
        {
          t.remove(t.vertices_begin());
          CGAL_assertion(t.is_valid(true));
        }
    }

  Random random(1284141159);
  std::cout << "Seed: " << random.get_seed () << std::endl;
  Random_points_in_square g(0.495, random);

  CGAL_assertion(t.is_valid());

  std::cout << "Removing first point" << std::endl;
  Vertex_handle vh0 = t.insert(Point(0.5, 0.5));
  CGAL_assertion(t.is_valid());
  t.remove(vh0);
  CGAL_assertion(t.is_valid());
  CGAL_assertion(t.empty());

  {
    Random random(1284141159);
    std::cout << "Seed: " << random.get_seed () << std::endl;

    Random_points_in_square g(0.495, random);
    Vector midpoint(0.5, 0.5);

    Triangulation t;
    CGAL_assertion(t.is_valid());

    std::cout << "Removing first point" << std::endl;
    Vertex_handle vh0 = t.insert(Point(0.5, 0.5));
    CGAL_assertion(t.is_valid());
    t.remove(vh0);
    CGAL_assertion(t.is_valid());
    CGAL_assertion(t.empty());

    std::cout << "Inserting random points and removing them." << std::endl;

    for (int i = 0; i < N_PTS; ++i)
      {
        t.insert(*(++g) + midpoint);
      }
    CGAL_assertion(t.is_valid());

    for (int i = 0; i < N_PTS; ++i)
      {
        // Find a random vertex
        Vertex_handle vh = t.locate(*(++g) + midpoint)->vertex(0);
        vh = t.get_original_vertex(vh);
        t.remove(vh);
        CGAL_assertion(t.is_valid());
      }
  }
  return 0;
}
