#include <CGAL/Simple_cartesian.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>

// Just to try and create ambiguities
#include <CGAL/boost/graph/io.h>
#include <CGAL/IO/io.h>

#include <CGAL/property_map.h>

#include <vector>
#include <cassert>
#include <string>
#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;

const double epsilon = 5e-4;

template < typename Type >
bool approx_equal_nt(const Type &t1, const Type &t2)
{
  if(t1 == t2)
    return true;
  if(CGAL::abs(t1 - t2) / (CGAL::max)(CGAL::abs(t1), CGAL::abs(t2)) < std::abs(t1)*epsilon)
    return true;

  std::cout << " Approximate comparison failed between : " << t1 << "  and  " << t2 << std::endl;
  std::cout << "abs(t1 - t2) = " <<CGAL::abs(t1 - t2) << "n";
  std::cout << "abs(t1) = " <<CGAL::abs(t1) << "n";
  std::cout << "abs(t2) = " <<CGAL::abs(t2) << "n";
  std::cout << "std::abs(t1)*epsilon = " <<std::abs(t1) * epsilon << "n";
  std::cout << "CGAL::abs(t1 - t2) / (CGAL::max)(CGAL::abs(t1), CGAL::abs(t2) == "
            << CGAL::abs(t1 - t2) / (CGAL::max)(CGAL::abs(t1), CGAL::abs(t2)) << std::endl;
  return false;
}

typedef std::pair<Point_3, Vector_3> PointVectorPair;
bool ps_are_equal(const CGAL::Point_set_3<Point_3, Vector_3>& ps,
                  const CGAL::Point_set_3<Point_3, Vector_3>& ps2)
{
  if(ps.size() != ps2.size())
    return false;

  typedef CGAL::Point_set_3<Point_3, Vector_3>::const_iterator Iterator;
  for(Iterator it1 = ps.begin(), it2 = ps2.begin(); it1!=ps.end() && it2!=ps2.end(); ++it1, ++it2)
  {
    const Point_3& p1 = ps.point(*it1);
    const Point_3& p2 = ps2.point(*it2);
    if(!(approx_equal_nt(p1.x(), p2.x()) && approx_equal_nt(p1.y(), p2.y()) && approx_equal_nt(p1.z(), p2.z())))
      return false;
  }

  return true;
}

typedef std::pair<Point_3, Vector_3> PointVectorPair;
bool points_are_equal(const std::vector<Point_3>& ps,
                      const std::vector<Point_3>& ps2)
{
  if(ps.size() != ps2.size())
    return false;

  typedef std::vector<Point_3>::const_iterator Iterator;
  for(Iterator it1 = ps.begin(), it2 = ps2.begin(); it1!=ps.end() && it2!=ps2.end(); ++it1, ++it2)
  {
    const Point_3& p1 = *it1;
    const Point_3& p2 = *it2;
    if(!(approx_equal_nt(p1.x(), p2.x()) && approx_equal_nt(p1.y(), p2.y()) && approx_equal_nt(p1.z(), p2.z())))
      return false;
  }

  return true;
}

bool test_points_with_np(std::string s)
{
  std::vector<PointVectorPair> points;
  points.push_back(std::make_pair(Point_3(0,0,0), Vector_3(0,0,0)));
  points.push_back(std::make_pair(Point_3(1,0,0), Vector_3(1,0,0)));
  points.push_back(std::make_pair(Point_3(0,1,0), Vector_3(0,1,0)));
  points.push_back(std::make_pair(Point_3(0,0,1), Vector_3(0,0,1)));
  points.push_back(std::make_pair(Point_3(1,1,1), Vector_3(1,1,1)));

  bool ok = CGAL::IO::write_points(s, points,
                                   CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                                                    .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
  assert(ok);
  std::vector<PointVectorPair> pv_pairs;
  ok = CGAL::IO::read_points(s, std::back_inserter(pv_pairs),
                             CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                                              .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
  assert(ok);
  assert(pv_pairs[0] == std::make_pair(Point_3(0,0,0), Vector_3(0,0,0)));
  assert(pv_pairs[1] == std::make_pair(Point_3(1,0,0), Vector_3(1,0,0)));
  assert(pv_pairs[2] == std::make_pair(Point_3(0,1,0), Vector_3(0,1,0)));
  assert(pv_pairs[3] == std::make_pair(Point_3(0,0,1), Vector_3(0,0,1)));
  assert(pv_pairs[4] == std::make_pair(Point_3(1,1,1), Vector_3(1,1,1)));

  return true;
}

#define CGAL_DEF_TEST_POINT_SET_3_FUNCTION(TYPE, type)                                             \
void test_##TYPE(const std::string& s)                                                             \
{                                                                                                  \
  std::cout << "Test Point_set_3: " << s << " extension: " << #TYPE << std::endl;                  \
  CGAL::Point_set_3<Point_3, Vector_3> ps;                                                         \
  bool ok = CGAL::IO::read_##TYPE(s, ps);                                                              \
  assert(ok);                                                                                      \
  ps.clear();                                                                                      \
  ok = CGAL::IO::read_##TYPE(s.c_str(), ps);                                                           \
  assert(ok);                                                                                      \
  ps.clear();                                                                                      \
  std::ifstream in(s);                                                                             \
  ok = CGAL::IO::read_##TYPE(in, ps);                                                                  \
  assert(ok);                                                                                      \
  const char* ext = type;                                                                          \
  std::string fname = "tmp.";                                                                      \
  fname.append(ext);                                                                               \
  ok = CGAL::IO::write_##TYPE(fname, ps);                                                              \
  assert(ok);                                                                                      \
  ok = CGAL::IO::write_##TYPE(fname, ps, CGAL::parameters::stream_precision(10));                      \
  assert(ok);                                                                                      \
  ok = CGAL::IO::write_##TYPE(fname.c_str(), ps);                                                      \
  assert(ok);                                                                                      \
  ok = CGAL::IO::write_##TYPE(fname.c_str(), ps, CGAL::parameters::stream_precision(10));              \
  assert(ok);                                                                                      \
  std::ofstream out(fname);                                                                        \
  ok = CGAL::IO::write_##TYPE(out, ps);                                                                \
  assert(ok);                                                                                      \
  std::ofstream out2(fname);                                                                       \
  ok = CGAL::IO::write_##TYPE(out2, ps, CGAL::parameters::stream_precision(10));                       \
  assert(ok);                                                                                      \
  CGAL::Point_set_3<Point_3, Vector_3> ps2;                                                        \
  std::ifstream is(fname);                                                                         \
  ok = CGAL::IO::read_##TYPE(is, ps2);                                                                 \
  assert(ok);                                                                                      \
  assert(ps_are_equal(ps, ps2));                                                                   \
  ok = CGAL::IO::write_point_set(fname, ps2);                                                          \
  assert(ok);                                                                                      \
  ok = CGAL::IO::write_point_set(fname, ps2, CGAL::parameters::stream_precision(10));                  \
  assert(ok);                                                                                      \
  ok = CGAL::IO::write_point_set(fname.c_str(), ps2);                                                  \
  assert(ok);                                                                                      \
  ok = CGAL::IO::write_point_set(fname.c_str(), ps2, CGAL::parameters::stream_precision(10));          \
  assert(ok);                                                                                      \
  ps2.clear();                                                                                     \
  ok = CGAL::IO::read_point_set(fname, ps2);                                                           \
  assert(ok);                                                                                      \
  assert(ps_are_equal(ps, ps2));                                                                   \
  }

CGAL_DEF_TEST_POINT_SET_3_FUNCTION(XYZ, "xyz")
CGAL_DEF_TEST_POINT_SET_3_FUNCTION(OFF, "off")
CGAL_DEF_TEST_POINT_SET_3_FUNCTION(PLY, "ply")
#ifdef CGAL_LINKED_WITH_LASLIB
void test_LAS(const std::string& s)
{
  std::cout << "Test Point_set_3: " << s << " extension: las" <<std::endl;
  CGAL::Point_set_3<Point_3, Vector_3> ps;
  bool ok = CGAL::IO::read_LAS(s, ps);
  assert(ok);
  ps.clear();
  ok = CGAL::IO::read_LAS(s.c_str(), ps);
  assert(ok);
  ps.clear();
  std::ifstream in(s, std::ios::binary);
  CGAL::IO::set_mode(in, CGAL::IO::BINARY);
  ok = CGAL::IO::read_LAS(in, ps);
  assert(ok);
  const char* ext = "las";
  std::string fname = "tmp.";
  fname.append(ext);
  ok = CGAL::IO::write_LAS(fname, ps);
  assert(ok);
  ok = CGAL::IO::write_LAS(fname.c_str(), ps);
  assert(ok);
  std::ofstream out(fname, std::ios::binary);
  CGAL::IO::set_mode(out, CGAL::IO::BINARY);
  ok = CGAL::IO::write_LAS(out, ps);
  assert(ok);
  CGAL::Point_set_3<Point_3, Vector_3> ps2;

  std::ifstream is(fname, std::ios::binary);
  CGAL::IO::set_mode(is, CGAL::IO::BINARY);
  ok = CGAL::IO::read_LAS(is, ps2);
  assert(ok);
  assert(ps_are_equal(ps, ps2));
  ok = CGAL::IO::write_point_set(fname, ps2);
  assert(ok);
  ok = CGAL::IO::write_point_set(fname, ps2, CGAL::parameters::stream_precision(10));
  assert(ok);
  ok = CGAL::IO::write_point_set(fname.c_str(), ps2);
  assert(ok);
  ok = CGAL::IO::write_point_set(fname.c_str(), ps2, CGAL::parameters::stream_precision(10));
  assert(ok);
  ps2.clear();
  ok = CGAL::IO::read_point_set(fname, ps2);
  assert(ok);
  assert(ps_are_equal(ps, ps2));
  }
#endif

#undef CGAL_DEF_INITIALIZE_ID_FUCNTION

#define CGAL_DEF_TEST_POINTS_FUNCTION(TYPE, type)                                                  \
void test_points_##TYPE(const std::string& s)                                                      \
{                                                                                                  \
  std::cout << "Test points: " << s << " extension: " << #TYPE << std::endl;                       \
  std::vector<Point_3> ps;                                                                         \
  bool ok = CGAL::IO::read_##TYPE(s, std::back_inserter(ps));                                          \
  assert(ok);                                                                                      \
  ps.clear();                                                                                      \
  ok = CGAL::IO::read_##TYPE(s.c_str(), std::back_inserter(ps));                                       \
  assert(ok);                                                                                      \
  ps.clear();                                                                                      \
  std::ifstream in(s);                                                                             \
  ok = CGAL::IO::read_##TYPE(in, std::back_inserter(ps));                                              \
  assert(ok);                                                                                      \
  const char* ext = type;                                                                          \
  std::string fname = "tmp.";                                                                      \
  fname.append(ext);                                                                               \
  ok = CGAL::IO::write_##TYPE(fname, ps);                                                              \
  assert(ok);                                                                                      \
  ok = CGAL::IO::write_##TYPE(fname, ps, CGAL::parameters::stream_precision(10));                      \
  assert(ok);                                                                                      \
  ok = CGAL::IO::write_##TYPE(fname.c_str(), ps);                                                      \
  assert(ok);                                                                                      \
  ok = CGAL::IO::write_##TYPE(fname.c_str(), ps, CGAL::parameters::stream_precision(10));              \
  assert(ok);                                                                                      \
  std::ofstream out(fname);                                                                        \
  ok = CGAL::IO::write_##TYPE(out, ps);                                                                \
  assert(ok);                                                                                      \
  std::ofstream out2(fname);                                                                       \
  ok = CGAL::IO::write_##TYPE(out2, ps, CGAL::parameters::stream_precision(10));                       \
  assert(ok);                                                                                      \
  std::vector<Point_3> ps2;                                                                        \
  std::ifstream is(fname);                                                                         \
  ok = CGAL::IO::read_##TYPE(is, std::back_inserter(ps2));                                             \
  assert(ok);                                                                                      \
  assert(points_are_equal(ps, ps2));                                                               \
  ok = CGAL::IO::write_points(fname, ps2);                                                             \
  assert(ok);                                                                                      \
  ps2.clear();                                                                                     \
  ok = CGAL::IO::read_points(fname, std::back_inserter(ps2));                                          \
  assert(ok);                                                                                      \
  assert(points_are_equal(ps, ps2));                                                               \
  }

CGAL_DEF_TEST_POINTS_FUNCTION(XYZ, "xyz")
CGAL_DEF_TEST_POINTS_FUNCTION(OFF, "off")
CGAL_DEF_TEST_POINTS_FUNCTION(PLY, "ply")
#ifdef CGAL_LINKED_WITH_LASLIB
void test_points_LAS(const std::string& s)
{
  std::cout << "Test points: " << s << " extension: LAS "<< std::endl;
  std::vector<Point_3> ps;
  bool ok = CGAL::IO::read_LAS(s, std::back_inserter(ps));
  assert(ok);
  ps.clear();
  ok = CGAL::IO::read_LAS(s.c_str(), std::back_inserter(ps));
  assert(ok);
  ps.clear();
  std::ifstream in(s, std::ios::binary);
  ok = CGAL::IO::read_LAS(in, std::back_inserter(ps));
  assert(ok);
  const char* ext = "las";
  std::string fname = "tmp.";
  fname.append(ext);
  ok = CGAL::IO::write_LAS(fname, ps);
  assert(ok);
  ok = CGAL::IO::write_LAS(fname, ps, CGAL::parameters::stream_precision(10));
  assert(ok);
  ok = CGAL::IO::write_LAS(fname.c_str(), ps);
  assert(ok);
  ok = CGAL::IO::write_LAS(fname.c_str(), ps, CGAL::parameters::stream_precision(10));
  assert(ok);
  std::ofstream out(fname, std::ios::binary);
  ok = CGAL::IO::write_LAS(out, ps);
  assert(ok);
  std::ofstream out2(fname, std::ios::binary);
  ok = CGAL::IO::write_LAS(out2, ps, CGAL::parameters::stream_precision(10));
  assert(ok);
  std::vector<Point_3> ps2;
  std::ifstream is(fname, std::ios::binary);
  ok = CGAL::IO::read_LAS(is, std::back_inserter(ps2));
  assert(ok);
  assert(points_are_equal(ps, ps2));
  ok = CGAL::IO::write_points(fname, ps2);
  assert(ok);
  ps2.clear();
  ok = CGAL::IO::read_points(fname, std::back_inserter(ps2));
  assert(ok);
  assert(points_are_equal(ps, ps2));
  }
#endif

#undef CGAL_DEF_INITIALIZE_ID_FUCNTION

int main()
{
  test_XYZ("data/read_test/ok_2.xyz");
  test_OFF("data/read_test/ok_1.off");
  test_PLY("data/read_test/simple_ascii.ply");
  test_PLY("data/read_test/simple.ply");

  test_points_XYZ("data/read_test/ok_2.xyz");
  test_points_OFF("data/read_test/ok_1.off");
  test_points_PLY("data/read_test/simple_ascii.ply");
  test_points_PLY("data/read_test/simple.ply");

#ifdef CGAL_LINKED_WITH_LASLIB
  test_LAS("data/read_test/pig_points.las");
  test_points_LAS("data/read_test/pig_points.las");
#endif

  test_points_with_np("test.xyz");
  test_points_with_np("test.off");
  test_points_with_np("test.ply");
  // don't test normals with LAS as it is not supported by the format

  return EXIT_SUCCESS;
}
