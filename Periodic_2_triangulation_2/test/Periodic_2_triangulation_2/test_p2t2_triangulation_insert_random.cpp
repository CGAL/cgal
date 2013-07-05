// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include "./types.h"

#ifndef N_PTS
#define N_PTS 10000
#endif

int main()
{
  Triangulation t;

  Random random(1284141159);
  std::cout << "Seed: " << random.get_seed () << std::endl;

  // Random_points_in_square g(0.495, random);
  Random_points_on_circle g(0.495, random);
  Vector midpoint(0.5, 0.5);

  for (int i = 0; i < N_PTS; ++i)
    {
      t.insert(*(++g) + midpoint);
    }

  if (!t.is_valid(true))
    {
      std::cout << "l:" << __LINE__ << std::endl;
      std::exit(1);
    }

  return 0;
}
