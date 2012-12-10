#include "./types.h"

#ifndef N_PTS
#define N_PTS 10000
#endif

int main() {
  Triangulation t;

  Random random(1284141159);
  std::cout << "Seed: " << random.get_seed () << std::endl;

  Random_points_in_square g(0.495, random);
  Vector midpoint(0.5, 0.5);
  
  for (int i = 0; i < N_PTS; ++i) {
    if (i%1000 == 0)
    {
      std::cout << "i: "<< i << std::endl;
      CGAL_assertion(t.is_valid());
    }
    t.insert(*(++g) + midpoint);
  }
  CGAL_assertion(t.is_valid());

  return 0;
}
