#include <limits>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/hierarchy_simplify_point_set.h>

#include <vector>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;

void test (std::vector<Point>& input,
	   std::ptrdiff_t result0 = 1, int result1 = 1, int result2 = 1, int result3 = 1, int result4 = 1)
{
  std::vector<Point>::iterator it = 
    CGAL::hierarchy_simplify_point_set (input.begin (), input.end (), 1);
  if (result0 > 0 && std::distance (input.begin (), it) != result0)
    exit (EXIT_FAILURE);

  it = CGAL::hierarchy_simplify_point_set (input.begin (), input.end ());
  if (result1 > 0 && std::distance (input.begin (), it) != result1)
    exit (EXIT_FAILURE);

  it = CGAL::hierarchy_simplify_point_set (input.begin (), input.end (), 100);
  if (result2 > 0 && std::distance (input.begin (), it) != result2)
    exit (EXIT_FAILURE);


  it = CGAL::hierarchy_simplify_point_set (input.begin (), input.end (), 1000, 0.1);
  if (result3 > 0 && std::distance (input.begin (), it) != result3)
    exit (EXIT_FAILURE);


  it = CGAL::hierarchy_simplify_point_set (input.begin (), input.end (),
					   CGAL::Identity_property_map<Point>(),
					   (std::numeric_limits<unsigned int>::max)(),
					   0.0001);
  if (result4 > 0 && std::distance (input.begin (), it) != result4)
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
  test (input, 2);

  // Test line
  for (std::size_t i = 0; i < 1000; ++ i)
    input.push_back (Point (0., 0., (double)i));
  test (input, input.size (), 128, 16, 1, 1);
  
  // Test plane
  for (std::size_t i = 0; i < 128; ++ i)
    for (std::size_t j = 0; j < 128; ++ j)
      input.push_back (Point (0., (double)j, (double)i));
  test (input, input.size (), 2048, 256, 32, 1);
  
  // Test random
  for (std::size_t i = 0; i < 10000; ++ i)
    input.push_back (Point (rand() / (FT)RAND_MAX,
  			    rand() / (FT)RAND_MAX,
  			    rand() / (FT)RAND_MAX));
  test (input, input.size (), -1, -1, -1, -1);
  
  return EXIT_SUCCESS;
}

