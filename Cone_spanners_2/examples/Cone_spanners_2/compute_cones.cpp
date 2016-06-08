#include <cstdlib>
#include <iostream>
#include <iterator>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_root_of.h>
#include <CGAL/Compute_cone_boundaries_2.h>

// select the kernel type
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_root_of   Kernel;
typedef Kernel::Point_2                   Point_2;
typedef Kernel::Direction_2               Direction_2;

int main(int argc, char ** argv) {

  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <no. of cones> [<direction-x> <direction-y>]" << std::endl;
    return 1;
  }

  unsigned int k = atoi(argv[1]);
  if (k<2) {
    std::cout << "The number of cones should be larger than 1!" << std::endl;
    return 1;
  }

  Direction_2 initial_direction;
  if (argc == 2)
    initial_direction = Direction_2(1, 0);  // default initial_direction
  else if (argc == 4)
    initial_direction = Direction_2(atof(argv[2]), atof(argv[3]));
  else {
    std::cout << "Usage: " << argv[0] << " <no. of cones> [<direction-x> <direction-y>]" << std::endl;
    return 1;
  }

  // construct the functor
  CGAL::Compute_cone_boundaries_2<Kernel> cones;
  // create the vector rays to store the results
  std::vector<Direction_2> rays(k);
  // compute the cone boundaries and store them in rays
  cones(k, initial_direction, rays.begin());

  // display the computed rays, starting from the initial direction, ccw order
  for (unsigned int i=0; i<k; i++)
    std::cout << "Ray " << i << ": " << rays[i] << std::endl;

  return 0;
}
