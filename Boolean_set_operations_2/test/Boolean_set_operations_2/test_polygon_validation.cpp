/* test file for polygon validation. Intended for testing the global functions
 * defined at Gps_polygon_validation.h
 */

#include <CGAL/assertions_behaviour.h>

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <iterator>
#include <string>
#include <sstream>
#include <iostream>

typedef CGAL::Exact_rational                       Number_type;
typedef CGAL::Cartesian<Number_type>               Kernel;
typedef CGAL::Gps_segment_traits_2<Kernel>         Traits_2;
typedef Traits_2::Polygon_2                        Polygon_2;
typedef Traits_2::Polygon_with_holes_2             Polygon_with_holes_2;

/* test files:
 *  1.  test1.dat---invalid polygon with holes. The hole is relatively simple
 *                  instead of strictly simple.
 *  2.  test2.dat---invalid polygon with holes. The hole overlaps the outer
 *                  boundary (the intersection results in a polygon).
 *  3.  test3.dat---invalid polygon with holes. Two holes intersect (the
 *                  intersection results in a polygon).
 *  4.  test4.dat---invalid polygon with holes. Two holes intersect (one
 *                  contains the other).
 *  5.  test5.dat---invalid polygon with holes. Two holes share an edge. (non
 *                  regularized intersection results in an edge).
 *  6.  test6.dat---invalid polygon with holes. A hole and the outer boundary
 *                  share an edge. (non regularized intersection results in an
 *                  edge).
 *  7.  test7.dat---invalid polygon with holes. The outer boundary is not
 *                  relatively simple because a "crossover" occurs at an
 *                  intersection
 *  8.  test8.dat---valid polygon with holes. Outer boundary is relatively
 *                  simple, case 1.
 *  9.  test9.dat---valid polygon with holes. Outer Boundary and holes are
 *                  pairwise disjoint except on vertices
 * 10. test10.dat---valid polygon with holes. Outer boundary is relatively
 *                  simple, case 2.
 *
 * test an input file. is_valid indicates the input polygon is valid. ErrorMsg
 * is displayed if the validation result does't equal is_valid.
 */

bool test(const std::string& filename, std::ofstream& outfile)
{
  std::ifstream input_file(filename.c_str());

  if (! input_file.is_open()) {
    std::cerr << "Failed to open the " << filename << std::endl;
    return (false);
  }
  // Read a polygon with holes from a file.
  Polygon_2 outer_pgn;
  size_t num_holes;
  bool is_valid;

  input_file >> outer_pgn;
  input_file >> num_holes;
  std::vector<Polygon_2> holes(num_holes);
  for (size_t k = 0; k < num_holes; ++k) input_file >> holes[k];
  input_file >> is_valid;

  Polygon_with_holes_2 P(outer_pgn, holes.begin(), holes.end());
  Traits_2 tr;
  bool test_valid = CGAL::is_valid_polygon_with_holes(P, tr);
  bool res = true;
  if (test_valid != is_valid) {
    res = false;
    outfile << "Error validating " << filename <<std::endl;
    //outfile << "P = " ;
    //print_polygon_with_holes (P);
    outfile << std::endl;
  }
  input_file.close();
  return res;
}

void special_warnings(const char*, const char* expr, const char* file, int line,
                      const char* msg)
{
  std::cerr << "  // CGAL: check violation! THIS MESSAGE IS PROBABLY WANTED."
            << std::endl
            << "  // Expression : " << expr << std::endl
            << "  // File       : " << file << std::endl
            << "  // Line       : " << line << std::endl
            << "  // Explanation: " << msg << std::endl
            << "  // Refer to the bug-reporting instructions at https://www.cgal.org/bug_report.html"
            << std::endl;
}

int main()
{
  std::cerr << "Modify the w-a-r-n-i-n-g-s handler...\n";
  CGAL::set_warning_handler(special_warnings);
  std::string testfile_prefix("data/validation/test");
  std::string testfile_suffix(".dat");
  const char* output_filename("data/validation/validation_test_output.txt");
  std::ofstream output_file (output_filename);
  if (! output_file.is_open()) {
    std::cerr << "Failed to open the " << output_filename <<std::endl;
    return (1);
  }

  int result = 0;
  for (size_t i = 1; i < 11; ++i) {
    std::stringstream strs;
    std::string si;
    strs << i;
    strs >> si;
    std::string filename = testfile_prefix + si + testfile_suffix;
    bool res =  test(filename, output_file);
    if (!res) {
      std::cout << "test " << i << " failed" << std::endl;
      result = 1;
    }
    else {
      std::cout <<"test " << i << " succeeded" << std::endl;
    }
  }
  if (result != 0) {
    std::cout << "SOME TESTS FAILED" << std::endl;
    return 1;
  }
  std::cout << "ALL TESTS SUCCEEDED!" << std::endl;
  return 0;
}
