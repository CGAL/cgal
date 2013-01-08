#include "./types.h"

#ifndef N_PTS
#define N_PTS 10000
#endif

int main() {
  Triangulation t;

  Random random(1284141159);
  std::cout << "Seed: " << random.get_seed () << std::endl;

  // Random_points_in_square g(0.495, random);
  Random_points_on_circle g(0.495, random);
  Vector midpoint(0.5, 0.5);

  int iterations = 0;
  while (1) {
    for (int i = 0; i < N_PTS; ++i) {
      t.insert(*(++g) + midpoint);
    }
      
    std::vector<Vertex_handle> vhs;
    for (Triangulation::Unique_vertex_iterator it = t.unique_vertices_begin();
         it != t.unique_vertices_end(); ++it) {
      vhs.push_back(it);
    }
      
    std::random_shuffle(vhs.begin(), vhs.end());
    if (iterations%2 == 0) {
      vhs.resize(vhs.size()/2);
    }
      
    for (size_t i=0; i<vhs.size(); ++i) {
      t.remove(vhs[i]); 
    }

    if (!t.is_valid(true)) {
      std::cout << "l:" << __LINE__ << std::endl;
      std::exit(1);
    }
    std::cout << "iteration: " << iterations++ << std::endl;
  }
  return 0;
}
