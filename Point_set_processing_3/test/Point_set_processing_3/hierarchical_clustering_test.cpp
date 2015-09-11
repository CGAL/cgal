#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/hierarchical_clustering.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;

void test (std::vector<Point>& input,
	   int result1 = 1, int result2 = 1, int result3 = 1, int result4 = 1)
{
  std::vector<Point> output;
  
  CGAL::hierarchical_clustering (input.begin (), input.end (),
				 std::back_inserter (output));
  if (result1 > 0 && output.size () != static_cast<std::size_t>(result1))
    exit (EXIT_FAILURE);
  output.clear ();
    
  CGAL::hierarchical_clustering (input.begin (), input.end (),
				 std::back_inserter (output), 100);
  if (result2 > 0 && output.size () != static_cast<std::size_t>(result2))
    exit (EXIT_FAILURE);
  output.clear ();
  
  CGAL::hierarchical_clustering (input.begin (), input.end (),
				 std::back_inserter (output), 1000, 0.1);
  if (result3 > 0 && output.size () != static_cast<std::size_t>(result3))
    exit (EXIT_FAILURE);
  output.clear ();

  CGAL::hierarchical_clustering (input.begin (), input.end (),
				 CGAL::Identity_property_map<Point>(),
				 std::back_inserter (output),
				 std::numeric_limits<unsigned int>::max(),
				 0.0001);
  if (result4 > 0 && output.size () != static_cast<std::size_t>(result4))
    exit (EXIT_FAILURE);

  input.clear ();
}


int main(void)
{

  std::vector<Point> input;

  // Test 1 point
  input.push_back (Point (0., 0., 0.));
  test (input);

  // Test twice the same point
  input.push_back (Point (0., 0., 0.));
  input.push_back (Point (0., 0., 0.));
  test (input);

  // Test 2 points
  input.push_back (Point (0., 0., 0.));
  input.push_back (Point (1., 0., 0.));
  test (input);

  // Test line
  for (std::size_t i = 0; i < 1000; ++ i)
    input.push_back (Point (0., 0., i));
  test (input, 128, 16, 1, 1);
  
  // Test plane
  for (std::size_t i = 0; i < 128; ++ i)
    for (std::size_t j = 0; j < 128; ++ j)
      input.push_back (Point (0., j, i));
  test (input, 2048, 256, 32, 1);
  
  // Test random
  for (std::size_t i = 0; i < 10000; ++ i)
    input.push_back (Point (rand() / (FT)RAND_MAX,
			    rand() / (FT)RAND_MAX,
			    rand() / (FT)RAND_MAX));
  test (input, -1, -1, -1, -1);
  
  return EXIT_SUCCESS;
}

