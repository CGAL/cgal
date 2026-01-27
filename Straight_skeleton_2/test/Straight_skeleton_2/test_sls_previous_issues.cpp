#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>
#include <CGAL/draw_straight_skeleton_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                   Point ;
typedef CGAL::Polygon_2<K>           Polygon_2 ;
typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef std::shared_ptr<Ss> SsPtr ;

// Issue #39 on CB
Polygon_2 data_0()
{
  Polygon_2 poly ;

  poly.push_back(Point(320, 1120));
  poly.push_back(Point(-141.55499021951055738, 1658.4808219227620611));
  poly.push_back(Point(-1115.9820207544369168, 1835.6551586833445526));
  poly.push_back(Point(760, 520));
  poly.push_back(Point(800, 560));
  poly.push_back(Point(840, 600));
  poly.push_back(Point(800,640));
  poly.push_back(Point(400, 1040));

  return poly ;
}

// Fixed by d62f376d6377cbc38120c82fd7db01acc0b3511d
Polygon_2 data_1()
{
  Polygon_2 poly ;

  poly.push_back(Point(355.7602676353456, 289.80005279541029));
  poly.push_back(Point(355.79700000000014, 289.47000000000014));
  poly.push_back(Point(355.79625000000016, 289.47075000000012));
  poly.push_back(Point(355.79587500000014, 289.47112500000014));
  poly.push_back(Point(355.79550000000017, 289.47150000000011));
  poly.push_back(Point(355.79512500000016, 289.47187500000013));
  poly.push_back(Point(355.79475000000014, 289.47225000000014));
  poly.push_back(Point(355.79437500000017, 289.47262500000011));
  poly.push_back(Point(355.79400000000015, 289.47300000000013));
  poly.push_back(Point(355.83765234375011, 289.15216015625015));
  poly.push_back(Point(355.88071875000014, 288.91628125000011));
  poly.push_back(Point(355.90653662109389, 288.82276806640635));
  poly.push_back(Point(355.93761328125015, 288.74157421875009));
  poly.push_back(Point(355.97575048828134, 288.66972607421883));
  poly.push_back(Point(356.02275000000014, 288.60425000000009));
  poly.push_back(Point(356.08041357421888, 288.5421723632814));
  poly.push_back(Point(356.15054296875013, 288.48051953125014));
  poly.push_back(Point(356.23493994140637, 288.4163178710939));
  poly.push_back(Point(356.33540625000018, 288.34659375000018));
  poly.push_back(Point(356.59175390625012, 288.17868359375012));
  poly.push_back(Point(356.93400000000014, 287.95300000000015));
  poly.push_back(Point(357.73257669067397, 287.87076025390638));
  poly.push_back(Point(358.41760180664073, 287.84929687500011));

  return poly ;
}

// # 10368
Polygon_2 data_2()
{
  Polygon_2 poly ;

  poly.push_back(Point(25.059928281249999, 25.810653964843748));
  poly.push_back(Point(25.555723457031252, 23.587245742187498));
  poly.push_back(Point(25.840255781250004, 23.678031562499996));
  poly.push_back(Point(26.693815253906251, 23.950390898437497));
  poly.push_back(Point(27.262852949218754, 24.131906757812498));

  return poly ;
}

void test(const Polygon_2& poly)
{
  SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());
  if(!iss)
  {
    std::cout << "Failed to create interior straight skeleton" << std::endl ;
    assert(false);
  }

  double lMaxOffset = 5 ;
  SsPtr oss = CGAL::create_exterior_straight_skeleton_2(lMaxOffset, poly);

  if(!oss)
  {
    std::cout << "Failed to create exterior straight skeleton" << std::endl ;
    assert(false);
  }

  CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
  CGAL::draw(*iss);

  CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*oss);
  CGAL::draw(*oss);
}

int main(int, char**)
{
  test(data_0());
  test(data_1());
  test(data_2());

  std::cout << "OK" << std::endl;

  return EXIT_SUCCESS ;
}
