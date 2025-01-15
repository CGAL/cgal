#define CGAL_SLS_TEST_SPEED_THINGS_UP_FOR_THE_TESTSUITE
#define CGAL_ENABLE_DISABLE_ASSERTIONS_AT_RUNTIME

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/create_weighted_offset_polygons_2.h>
#include <CGAL/create_weighted_offset_polygons_from_polygon_with_holes_2.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>

#include <boost/shared_ptr.hpp>

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel            EPECK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  EPECK_w_sqrt;

namespace CGAL {

template<typename K>
class Test_polygon_2 : public CGAL::Polygon_2<K> {
    typedef CGAL::Polygon_2<K> Base;
    Test_polygon_2(const Base&);
public:
    using Base::Base;
};

template<typename K>
class Test_polygon_with_holes_2 : public CGAL::Polygon_with_holes_2<K> {
    typedef CGAL::Polygon_with_holes_2<K> Base;
    Test_polygon_with_holes_2(const Base&);
public:
    using Base::Base;
};

} // namespace CGAL

using namespace CGAL;

template <typename K>
void test_API()
{
  typedef typename K::FT                                             FT;
  typedef typename K::Point_2                                        Point_2;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;

  typedef CGAL::Polygon_2<EPICK>                                     Polygon_2_EPICK;
  typedef CGAL::Polygon_with_holes_2<EPICK>                          Polygon_with_holes_2_EPICK;

  typedef CGAL::Test_polygon_2<K>                                    Test_Polygon_2;
  typedef CGAL::Test_polygon_with_holes_2<K>                         Test_Polygon_with_holes_2;

  typedef CGAL::Test_polygon_2<EPICK>                                Test_Polygon_2_EPICK;
  typedef CGAL::Test_polygon_with_holes_2<EPICK>                     Test_Polygon_with_holes_2_EPICK;

  std::vector<Point_2> v;
  Polygon_2 p;
  Polygon_with_holes_2 pwh;
  std::vector<std::vector<FT> > weights;

  std::vector< boost::shared_ptr<Polygon_2> > res;
  std::vector< boost::shared_ptr<Polygon_2_EPICK> > res_EPICK;
  std::vector< boost::shared_ptr<Polygon_with_holes_2> > res_wh;
  std::vector< boost::shared_ptr<Polygon_with_holes_2_EPICK> > res_wh_EPICK;

  std::vector< boost::shared_ptr<Test_Polygon_2> > res_test;
  std::vector< boost::shared_ptr<Test_Polygon_2_EPICK> > res_test_EPICK;
  std::vector< boost::shared_ptr<Test_Polygon_with_holes_2> > res_wh_test;
  std::vector< boost::shared_ptr<Test_Polygon_with_holes_2_EPICK> > res_wh_test_EPICK;

  // First kernel is the offset construction (and thus output kernel), second kernel is the skeleton construction

  // simple interior, no holes
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights) ;
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, EPICK()) ;
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, EPICK(), EPICK()) ;
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, EPICK(), K()) ;
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2_EPICK>(0.1, p, weights, EPICK(), EPICK()) ;
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2_EPICK>(0.1, p, weights, EPICK(), K()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, K()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, K(), EPICK()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, K(), K()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2(FT(0.1), p, weights, K(), K()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(0.1, p, weights, K(), EPICK()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(FT(0.1), p, weights, K(), EPICK()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(0.1, p, weights, K(), K()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(FT(0.1), p, weights, K(), K()) ;

  res_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, p, weights, EPICK(), EPICK()) ;
  res_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, p, weights, EPICK(), K()) ;
  res_test = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(0.1, p, weights, K(), K()) ;
  res_test = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(FT(0.1), p, weights, K(), K()) ;

  res_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, v, weights, EPICK(), EPICK()) ;
  res_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, v, weights, EPICK(), K()) ;
  res_test = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(0.1, v, weights, K(), K()) ;
  res_test = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(FT(0.1), v, weights, K(), K()) ;

  // simple interior, holes
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights) ;
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, EPICK()) ;
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, EPICK(), K()) ;
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2_EPICK>(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2_EPICK>(0.1, pwh, weights, EPICK(), K()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, K()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, K(), EPICK()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, K(), K()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2(FT(0.1), pwh, weights, K(), K()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(0.1, pwh, weights, K(), EPICK()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(FT(0.1), pwh, weights, K(), EPICK()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(0.1, pwh, weights, K(), K()) ;
  res = create_interior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(FT(0.1), pwh, weights, K(), K()) ;

  res_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, pwh, weights, EPICK(), K()) ;
  res_test = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(0.1, pwh, weights, K(), K()) ;
  res_test = create_interior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(FT(0.1), pwh, weights, K(), K()) ;

  // simple exterior, no holes
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights) ;
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, EPICK()) ;
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, EPICK(), EPICK()) ;
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, EPICK(), K()) ;
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2_EPICK>(0.1, p, weights, EPICK(), EPICK()) ;
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2_EPICK>(0.1, p, weights, EPICK(), K()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, K()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, K(), EPICK()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, p, weights, K(), K()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2(FT(0.1), p, weights, K(), K()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(0.1, p, weights, K(), EPICK()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(FT(0.1), p, weights, K(), EPICK()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(0.1, p, weights, K(), K()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(FT(0.1), p, weights, K(), K()) ;

  res_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, p, weights, EPICK(), EPICK()) ;
  res_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, p, weights, EPICK(), K()) ;
  res_test = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(0.1, p, weights, K(), K()) ;
  res_test = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(FT(0.1), p, weights, K(), K()) ;

  res_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, v, weights, EPICK(), EPICK()) ;
  res_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, v, weights, EPICK(), K()) ;
  res_test = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(0.1, v, weights, K(), K()) ;
  res_test = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(FT(0.1), v, weights, K(), K()) ;

  // simple exterior, holes
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights) ;
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, EPICK()) ;
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, EPICK(), K()) ;
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2_EPICK>(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2_EPICK>(0.1, pwh, weights, EPICK(), K()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, K()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, K(), EPICK()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2(0.1, pwh, weights, K(), K()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2(FT(0.1), pwh, weights, K(), K()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(0.1, pwh, weights, K(), EPICK()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(FT(0.1), pwh, weights, K(), EPICK()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(0.1, pwh, weights, K(), K()) ;
  res = create_exterior_weighted_skeleton_and_offset_polygons_2<Polygon_2>(FT(0.1), pwh, weights, K(), K()) ;

  res_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2_EPICK>(0.1, pwh, weights, EPICK(), K()) ;
  res_test = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(0.1, pwh, weights, K(), K()) ;
  res_test = create_exterior_weighted_skeleton_and_offset_polygons_2<Test_Polygon_2>(FT(0.1), pwh, weights, K(), K()) ;

  // Same, but the result has holes --------------------

  // arranged interior, no holes
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights) ;
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, EPICK()) ;
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, EPICK(), EPICK()) ;
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, EPICK(), K()) ;
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2_EPICK>(0.1, p, weights, EPICK(), EPICK()) ;
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2_EPICK>(0.1, p, weights, EPICK(), K()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, K()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, K(), EPICK()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, K(), K()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(FT(0.1), p, weights, K(), K()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(0.1, p, weights, K(), EPICK()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(FT(0.1), p, weights, K(), EPICK()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(0.1, p, weights, K(), K()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(FT(0.1), p, weights, K(), K()) ;

  res_wh_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, p, weights, EPICK(), EPICK()) ;
  res_wh_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, p, weights, EPICK(), K()) ;
  res_wh_test = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(0.1, p, weights, K(), K()) ;
  res_wh_test = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(FT(0.1), p, weights, K(), K()) ;

  res_wh_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, v, weights, EPICK(), EPICK()) ;
  res_wh_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, v, weights, EPICK(), K()) ;
  res_wh_test = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(0.1, v, weights, K(), K()) ;
  res_wh_test = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(FT(0.1), v, weights, K(), K()) ;

  // arranged interior, holes
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights) ;
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights, EPICK()) ;
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights, EPICK(), K()) ;
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2_EPICK>(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_wh_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2_EPICK>(0.1, pwh, weights, EPICK(), K()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights, K()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights, K(), EPICK()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(FT(0.1), pwh, weights, K(), K()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(0.1, pwh, weights, K(), EPICK()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(FT(0.1), pwh, weights, K(), EPICK()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(0.1, pwh, weights, K(), K()) ;
  res_wh = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(FT(0.1), pwh, weights, K(), K()) ;

  res_wh_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_wh_test_EPICK = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, pwh, weights, EPICK(), K()) ;
  res_wh_test = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(0.1, pwh, weights, K(), K()) ;
  res_wh_test = create_interior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(FT(0.1), pwh, weights, K(), K()) ;

  // arranged exterior, no holes
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights) ;
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, EPICK()) ;
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, EPICK(), EPICK()) ;
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, EPICK(), K()) ;
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2_EPICK>(0.1, p, weights, EPICK(), EPICK()) ;
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2_EPICK>(0.1, p, weights, EPICK(), K()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, K()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, K(), EPICK()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, p, weights, K(), K()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(0.1, p, weights, K(), EPICK()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(FT(0.1), p, weights, K(), EPICK()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(0.1, p, weights, K(), K()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(FT(0.1), p, weights, K(), K()) ;

  res_wh_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, p, weights, EPICK(), EPICK()) ;
  res_wh_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, p, weights, EPICK(), K()) ;
  res_wh_test = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(0.1, p, weights, K(), K()) ;
  res_wh_test = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(FT(0.1), p, weights, K(), K()) ;

  res_wh_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, v, weights, EPICK(), EPICK()) ;
  res_wh_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, v, weights, EPICK(), K()) ;
  res_wh_test = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(0.1, v, weights, K(), K()) ;
  res_wh_test = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(FT(0.1), v, weights, K(), K()) ;

  // arranged exterior, holes
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights) ;
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights, EPICK()) ;
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights, EPICK(), K()) ;
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2_EPICK>(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_wh_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2_EPICK>(0.1, pwh, weights, EPICK(), K()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights, K()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(0.1, pwh, weights, K(), EPICK()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(0.1, pwh, weights, K(), K()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(FT(0.1), pwh, weights, K(), EPICK()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(0.1, pwh, weights, K(), K()) ;
  res_wh = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Polygon_with_holes_2>(FT(0.1), pwh, weights, K(), K()) ;

  res_wh_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, pwh, weights, EPICK(), EPICK()) ;
  res_wh_test_EPICK = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2_EPICK>(0.1, pwh, weights, EPICK(), K()) ;
  res_wh_test = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(0.1, pwh, weights, K(), K()) ;
  res_wh_test = create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2<Test_Polygon_with_holes_2>(FT(0.1), pwh, weights, K(), K()) ;
}

template <typename K>
void test_kernel()
{
  void (*dummy_ptr)() = &test_API<K>;
}

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  test_kernel<EPICK>();
  test_kernel<EPECK>();
  test_kernel<EPECK_w_sqrt>();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
