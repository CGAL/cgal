// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include "./types.h"
#include <map>

void greedy_flip_long_edges(Triangulation &t)
{
  typedef std::multimap<double, std::pair<Vertex_handle, Vertex_handle> > Edge_map;

  Face_handle f;
  int i;

  Edge_map edge_map;
  for (Triangulation::Edge_iterator eit = t.edges_begin(); eit != t.edges_end(); ++eit)
    {
      double sqr_length_orig = t.segment(eit).squared_length();
      if (sqr_length_orig > 0.166)
        {
          f = eit->first;
          i = eit->second;
          edge_map.insert(std::make_pair(sqr_length_orig,
                                         std::make_pair(f->vertex(t.ccw(i)),
                                             f->vertex(t.cw(i)))));
        }
    }

  bool is_1_cover = t.is_1_cover();
  for (Edge_map::reverse_iterator it = edge_map.rbegin(); it != edge_map.rend(); ++it)
    {
      double sqr_length_orig = it->first;
      if (t.is_edge(it->second.first, it->second.second, f, i))
        {
          if (t.flippable(f, i))
            {
              t.flip(f, i);
              // We flipped enough long edges, when we go to the 1-cover all faces are invalidated
              if (is_1_cover != t.is_1_cover())
                return;
              double sqr_length_new = t.segment(f, t.ccw(i)).squared_length();
              if (sqr_length_orig < sqr_length_new)
                {
                  std::cout << sqr_length_orig << std::endl;
                  t.flip(f, t.ccw(i));
                }
            }
        }
    }
}

int main()
{
  Point p;
  Triangulation t;

  const int N = 4;

  // Insert the first point
  for (int y = 0; y < N; y++)
    {
      for (int x = 0; x < N; x++)
        {
          t.insert(Point((1.0 / N)*x, (1.0 / N)*y));
        }
    }

  greedy_flip_long_edges(t);
  greedy_flip_long_edges(t);
  CGAL_assertion(t.is_valid());

  size_t n_vertices = t.number_of_vertices();
  size_t n_faces = t.number_of_faces();
  CGAL_USE(n_vertices);
  CGAL_USE(n_faces);

  std::cout << "Converting to 1 cover" << std::endl;
  t.convert_to_1_sheeted_covering();
  std::cout << "... done" << std::endl;
  CGAL_assertion(t.is_1_cover());
  CGAL_assertion(t.is_valid());
  CGAL_assertion(t.number_of_vertices() == n_vertices);
  CGAL_assertion(t.number_of_faces() == n_faces);

  std::cout << "Converting to 9 cover" << std::endl;
  t.convert_to_9_sheeted_covering();
  std::cout << "... done" << std::endl;
  CGAL_assertion(!t.is_1_cover());
  CGAL_assertion(t.is_valid());
  CGAL_assertion(t.number_of_vertices() == n_vertices);
  CGAL_assertion(t.number_of_faces() == n_faces);

  std::cout << "Converting to 1 cover" << std::endl;
  t.convert_to_1_sheeted_covering();
  std::cout << "... done" << std::endl;
  CGAL_assertion(t.is_1_cover());
  CGAL_assertion(t.is_valid());
  CGAL_assertion(t.number_of_vertices() == n_vertices);
  CGAL_assertion(t.number_of_faces() == n_faces);

  return 0;
}
