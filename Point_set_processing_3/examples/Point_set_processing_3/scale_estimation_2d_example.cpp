#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/estimate_scale.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Random.h>

#include <vector>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_2 Point_2;

int main (int, char**)
{
  std::vector<Point_2> samples;
  samples.reserve (100000);

  // Generate circle with gradually variable noise
  //  - noise-free for points with x close to (-1)
  //  - noisy for points with x close to (+1)
  for (std::size_t i = 0; i < 100000; ++ i)
    {
      FT theta = CGAL::get_default_random().get_double(0, 2. * CGAL_PI);
      FT noise = 0.5 * (std::cos(theta) + 1.) * CGAL::get_default_random().get_double(0., 0.2);
      int mult = (CGAL::get_default_random().get_bool() ? 1 : -1);
      samples.push_back (Point_2 (std::cos(theta) * (1. + mult * noise * noise),
                                  std::sin(theta) * (1. + mult * noise * noise)));
    }

  // Search for local scales on 3 different locations
  std::vector<Point_2> queries;
  queries.reserve (3);
  queries.push_back (Point_2 (-1., 0.));
  queries.push_back (Point_2 (0., 1.));
  queries.push_back (Point_2 (1., 0.));

  std::vector<std::size_t> k_scales;
  CGAL::estimate_local_k_neighbor_scales (samples,
                                          queries,
                                          std::back_inserter (k_scales));

  // Display results
  std::cerr << "K-Scales found:" << std::endl
            << " - On noise-free region: " << k_scales[0] << std::endl // Should be small
            << " - On moderately noisy region: " << k_scales[1] << std::endl // Should be higher
            << " - On very noisy region: " << k_scales[2] << std::endl; // Should be even higher

  return EXIT_SUCCESS;
}

