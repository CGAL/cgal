#include <CGAL/Simple_cartesian.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/property_map.h>

#include <vector>
#include <deque>
#include <iostream>
#include <fstream>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef std::pair<Point_3, Vector_3> PointVectorPair;

// this is going to be custom OutputIterator value_type
struct dummy_counter {
  static std::size_t counter;

  dummy_counter() { ++counter; }
  operator std::size_t() { return counter-1; }
};
std::size_t dummy_counter::counter = 0;

bool check_points_and_vectors(
  const boost::vector_property_map<Point_3>& points,
  const boost::vector_property_map<Vector_3>& normals,
  const std::vector<PointVectorPair>& pv_pairs,
  const std::vector<std::size_t>& indices)
{
  if(pv_pairs.size() != indices.size()) {
    std::cerr << "Error: inconsistency between point / normal size." << std::endl;
    return false;
  }

  for(std::size_t i = 0; i < pv_pairs.size(); ++i ) {
    if(pv_pairs[i].first != points[i]) {
      std::cerr << "Error: points are not equal." << std::endl;
      return false;
    }
    if(pv_pairs[i].second != normals[i]) {
      std::cerr << "Error: normals are not equal." << std::endl;
      return false;
    }
  }
  return true;
}

bool check_points(
  const boost::vector_property_map<Point_3>& points_1,
  const std::vector<Point_3>& points_2,
  const std::vector<std::size_t>& indices) 
{
  if(points_2.size() != indices.size()) {
    std::cerr << "Error: inconsistency between point / normal size." << std::endl;
    return false;
  }

  for(std::size_t i = 0; i < points_2.size(); ++i ) {
    if(points_2[i] != points_1[i]) {
      std::cerr << "Error: points are not equal." << std::endl;
      return false;
    }
  }
  return true;
}

bool test_no_deduction_points_and_normals_xyz(const char* file_name)
{
  boost::vector_property_map<Point_3> points;
  boost::vector_property_map<Vector_3> normals;
  std::vector<std::size_t> indices;

  std::vector<PointVectorPair> pv_pairs;

  // read with custom output iterator type
  dummy_counter::counter = 0;
  std::ifstream input(file_name);
  CGAL::read_xyz_points_and_normals<dummy_counter>(
    input, back_inserter(indices), points, normals, Kernel());

  // read with ordinary pmaps
  input.clear();
  input.close();
  input.open(file_name);
  CGAL::read_xyz_points_and_normals(
    input, back_inserter(pv_pairs), 
    CGAL::First_of_pair_property_map<PointVectorPair>(), 
    CGAL::Second_of_pair_property_map<PointVectorPair>(),
    Kernel());

  return check_points_and_vectors(points, normals, pv_pairs, indices);
}

bool test_no_deduction_points_and_normals_off(const char* file_name)
{
  boost::vector_property_map<Point_3> points;
  boost::vector_property_map<Vector_3> normals;
  std::vector<std::size_t> indices;

  std::vector<PointVectorPair> pv_pairs;

  // read with custom output iterator type
  dummy_counter::counter = 0;
  std::ifstream input(file_name);
  CGAL::read_off_points_and_normals<dummy_counter>(
    input, back_inserter(indices), points, normals, Kernel());

  // read with ordinary pmaps
  input.clear();
  input.close();
  input.open(file_name);
  CGAL::read_off_points_and_normals(
    input, back_inserter(pv_pairs), 
    CGAL::First_of_pair_property_map<PointVectorPair>(), 
    CGAL::Second_of_pair_property_map<PointVectorPair>(),
    Kernel());

  return check_points_and_vectors(points, normals, pv_pairs, indices);
}

bool test_no_deduction_points_xyz(const char* file_name)
{
  boost::vector_property_map<Point_3> points_1; \
  std::vector<std::size_t> indices;

  std::vector<Point_3> points_2;

  // read with custom output iterator type
  dummy_counter::counter = 0;
  std::ifstream input(file_name);
  CGAL::read_xyz_points<dummy_counter>(
    input, back_inserter(indices), points_1, Kernel());

  // read with ordinary pmaps
  input.clear();
  input.close();
  input.open(file_name);
  CGAL::read_xyz_points(
    input, back_inserter(points_2), 
    CGAL::Identity_property_map<Point_3>(), 
    Kernel());

  return check_points(points_1, points_2, indices);
}

bool test_no_deduction_points_off(const char* file_name)
{
  boost::vector_property_map<Point_3> points_1;
  std::vector<std::size_t> indices;

  std::vector<Point_3> points_2;

  // read with custom output iterator type
  dummy_counter::counter = 0;
  std::ifstream input(file_name);
  CGAL::read_off_points<dummy_counter>(
    input, back_inserter(indices), points_1, Kernel());

  // read with ordinary pmaps
  input.clear();
  input.close();
  input.open(file_name);
  CGAL::read_off_points(
    input, back_inserter(points_2), 
    CGAL::Identity_property_map<Point_3>(), 
    Kernel());

  return check_points(points_1, points_2, indices);
}

void compile_test() {
  std::deque<Point_3> points;
  std::deque<Vector_3> normals;
  std::deque<PointVectorPair> pv_pairs;
  std::ifstream input;
  
  input.open("data/read_test/simple.xyz");
  CGAL::read_xyz_points(
    input,
    std::front_inserter(points));
  input.clear();
  input.close();
  
  input.open("data/read_test/simple.xyz");
  CGAL::read_xyz_points(
    input,
    std::front_inserter(points),
    CGAL::Identity_property_map<Point_3>());
  input.clear();
  input.close();

  input.open("data/read_test/simple.xyz");
  CGAL::read_xyz_points(
    input,
    std::front_inserter(points),
    CGAL::Identity_property_map<Point_3>(),
    Kernel());
  input.clear();
  input.close();

  // this will span all OutputIteratorValueType versions
  input.open("data/read_test/simple.xyz");
  CGAL::read_xyz_points<Point_3>(
    input,
    std::front_inserter(points));
  input.clear();
  input.close();
  //-----------------------------------------------------------------------
  input.open("data/read_test/simple.off");
  CGAL::read_off_points(
    input,
    std::front_inserter(points));
  input.clear();
  input.close();
  
  input.open("data/read_test/simple.off");
  CGAL::read_off_points(
    input,
    std::front_inserter(points),
    CGAL::Identity_property_map<Point_3>());
  input.clear();
  input.close();

  input.open("data/read_test/simple.off");
  CGAL::read_off_points(
    input,
    std::front_inserter(points),
    CGAL::Identity_property_map<Point_3>(),
    Kernel());
  input.clear();
  input.close();

  // this will span all OutputIteratorValueType versions
  input.open("data/read_test/simple.off");
  CGAL::read_off_points<Point_3>(
    input,
    std::front_inserter(points));
  input.clear();
  input.close();
  //-----------------------------------------------------------------------
  input.open("data/read_test/simple.xyz");
  CGAL::read_xyz_points_and_normals(
    input,
    std::front_inserter(points),
    boost::dummy_property_map());
  input.clear();
  input.close();

  input.open("data/read_test/simple.xyz");
  CGAL::read_xyz_points_and_normals(
    input,
    std::front_inserter(pv_pairs),
    CGAL::First_of_pair_property_map<PointVectorPair>(), 
    CGAL::Second_of_pair_property_map<PointVectorPair>());
  input.clear();
  input.close();

  input.open("data/read_test/simple.xyz");
  CGAL::read_xyz_points_and_normals(
    input,
    std::front_inserter(pv_pairs),
    CGAL::First_of_pair_property_map<PointVectorPair>(), 
    CGAL::Second_of_pair_property_map<PointVectorPair>(),
    Kernel());
  input.clear();
  input.close();

  input.open("data/read_test/simple.xyz");
  CGAL::read_xyz_points_and_normals<Point_3>(
    input,
    std::front_inserter(points),
    boost::dummy_property_map());
  input.clear();
  input.close();
  //-----------------------------------------------------------------------
  input.open("data/read_test/simple.off");
  CGAL::read_off_points_and_normals(
    input,
    std::front_inserter(points),
    boost::dummy_property_map());
  input.clear();
  input.close();

  input.open("data/read_test/simple.off");
  CGAL::read_off_points_and_normals(
    input,
    std::front_inserter(pv_pairs),
    CGAL::First_of_pair_property_map<PointVectorPair>(), 
    CGAL::Second_of_pair_property_map<PointVectorPair>());
  input.clear();
  input.close();

  input.open("data/read_test/simple.off");
  CGAL::read_off_points_and_normals(
    input,
    std::front_inserter(pv_pairs),
    CGAL::First_of_pair_property_map<PointVectorPair>(), 
    CGAL::Second_of_pair_property_map<PointVectorPair>(),
    Kernel());
  input.clear();
  input.close();

  input.open("data/read_test/simple.off");
  CGAL::read_off_points_and_normals<Point_3>(
    input,
    std::front_inserter(points),
    boost::dummy_property_map());
  input.clear();
  input.close();
}

int main() {
  if(!test_no_deduction_points_and_normals_xyz("data/read_test/simple.xyz")) {
    return EXIT_FAILURE;
  }
  std::cerr << "test_no_deduction_points_and_normals_xyz OK." << std::endl;

  if(!test_no_deduction_points_and_normals_off("data/read_test/simple.off")) {
    return EXIT_FAILURE;
  }
  std::cerr << "test_no_deduction_points_and_normals_off OK." << std::endl;

  if(!test_no_deduction_points_xyz("data/read_test/simple.xyz")) {
    return EXIT_FAILURE;
  }
  std::cerr << "test_no_deduction_points_xyz OK." << std::endl;

  if(!test_no_deduction_points_off("data/read_test/simple.off")) {
    return EXIT_FAILURE;
  }
  std::cerr << "test_no_deduction_points_off OK." << std::endl;

  compile_test();
  return EXIT_SUCCESS;
}
