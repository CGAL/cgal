#include "./types.h"

const int N_PTS = 50;

int main() {
  Triangulation t;

  Vector midpoint(0.5, 0.5);

  // Vertex_handle vh = t.nearest_vertex(Point(0,0) + midpoint);
  // CGAL_assertion(vh == Vertex_handle());

  Vertex_handle vh_midpoint = t.insert(Point(0,0) + midpoint);
  // vh = t.nearest_vertex(Point(0,0) + midpoint);
  // CGAL_assertion(vh == vh_midpoint);
  t.remove(vh_midpoint);
  CGAL_assertion(t.empty());

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
  
  std::cout << "Inserting random points and removing them." << std::endl;
  
  for (int i = 0; i < N_PTS; ++i) {
    t.insert(*(++g) + midpoint);
  }
  CGAL_assertion(t.is_valid());

  return 0;
  {
    // TODO(NGHK):
    CGAL_assertion(false && "Takes too long");
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

    for (int i = 0; i < N_PTS; ++i) {
      t.insert(*(++g) + midpoint);
    }
    CGAL_assertion(t.is_valid());

    for (int i = 0; i < N_PTS; ++i) {
      std::cout << "Remove: " << i << std::endl;
      // Find a random vertex
      //     Vertex_handle vh = t.locate(*(++g) + midpoint)->vertex(0);
      Vertex_handle vh = t.locate(Point(0.0, 0.0) + midpoint)->vertex(0);
      t.remove(vh);
      CGAL_assertion(t.is_valid());
    }
  }
  return 0;
}
