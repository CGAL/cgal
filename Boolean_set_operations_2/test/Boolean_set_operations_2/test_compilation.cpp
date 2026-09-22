
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_non_caching_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Polygon_set_2.h>

//typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;
typedef int Number_type;

typedef CGAL::Simple_cartesian<Number_type>             Kernel;

typedef CGAL::Gps_segment_traits_2<Kernel>              Traits;
typedef CGAL::Polygon_set_2<Kernel>                     Ps;

typedef CGAL::Arr_segment_traits_2<Kernel>              Arr_traits;
typedef CGAL::Gps_traits_2<Arr_traits>                  General_traits;
typedef CGAL::General_polygon_set_2<General_traits>     Gps;

typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Nc_traits;
typedef CGAL::Gps_segment_traits_2<Kernel,
                                   std::vector<Kernel::Point_2>,
                                   Nc_traits>           Traits_non_caching;
typedef CGAL::General_polygon_set_2<Traits_non_caching> Gps_non_caching;

template <class GPS>
void test()
{
  typedef typename GPS::Traits_2                        Traits;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::Polygon_2                    Polygon_2;
  typedef typename Traits::Polygon_with_holes_2         Polygon_with_holes_2;

  Polygon_2 pgn1, pgn2;
  Polygon_with_holes_2  pgn_with_holes1, pgn_with_holes2;
  std::vector<Polygon_2>             polygons;
  std::vector<Polygon_with_holes_2>  polygons_with_holes;
  GPS gps;
  GPS other;

  gps.intersection(pgn1);
  gps.intersection(pgn_with_holes1);
  gps.intersection(gps);

  gps.join(pgn1);
  gps.join(pgn_with_holes1);
  gps.join(gps);

  gps.difference(pgn1);
  gps.difference(pgn_with_holes1);
  gps.difference(gps);

  gps.symmetric_difference(pgn1);
  gps.symmetric_difference(pgn_with_holes1);
  gps.symmetric_difference(gps);

  gps.intersection(polygons.begin(), polygons.end());
  gps.intersection(polygons_with_holes.begin(), polygons_with_holes.end());
  gps.intersection(polygons.begin(), polygons.end(),
                   polygons_with_holes.begin(), polygons_with_holes.end());

  gps.join(polygons.begin(), polygons.end());
  gps.join(polygons_with_holes.begin(), polygons_with_holes.end());
  gps.join(polygons.begin(), polygons.end(),
           polygons_with_holes.begin(), polygons_with_holes.end());

  gps.symmetric_difference(polygons.begin(), polygons.end());
  gps.symmetric_difference(polygons_with_holes.begin(),
                           polygons_with_holes.end());
  gps.symmetric_difference(polygons.begin(), polygons.end(),
                           polygons_with_holes.begin(),
                           polygons_with_holes.end());

  gps.do_intersect(pgn1);
  gps.do_intersect(pgn_with_holes1);
  gps.do_intersect(other);

  gps.complement();
  gps.complement(gps);

  gps.do_intersect(polygons.begin(), polygons.end());
  gps.do_intersect(polygons_with_holes.begin(), polygons_with_holes.end());
  gps.do_intersect(polygons.begin(), polygons.end(),
                   polygons_with_holes.begin(), polygons_with_holes.end());

  std::vector<Polygon_with_holes_2>  result;

  Traits tr;
  // global functions

  CGAL::do_intersect(pgn1, pgn2);
  CGAL::do_intersect(pgn1, pgn_with_holes2);
  CGAL::do_intersect(pgn_with_holes1, pgn2);
  CGAL::do_intersect(pgn_with_holes1, pgn_with_holes2);
  CGAL::do_intersect(polygons.begin(), polygons.end());
  CGAL::do_intersect(polygons_with_holes.begin(), polygons_with_holes.end());
  CGAL::do_intersect(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end());


  CGAL::do_intersect(pgn1, pgn2, tr);
  CGAL::do_intersect(pgn1, pgn_with_holes2, tr);
  CGAL::do_intersect(pgn_with_holes1, pgn2, tr);
  CGAL::do_intersect(pgn_with_holes1, pgn_with_holes2, tr);
  CGAL::do_intersect(pgn_with_holes1, pgn_with_holes2);
  CGAL::do_intersect(polygons.begin(), polygons.end(), tr);
  CGAL::do_intersect(polygons_with_holes.begin(), polygons_with_holes.end(), tr);
  CGAL::do_intersect(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end(), tr);

  CGAL::intersection(pgn1, pgn2, std::back_inserter(result));
  CGAL::intersection(pgn1, pgn_with_holes2, std::back_inserter(result));
  CGAL::intersection(pgn_with_holes1, pgn2, std::back_inserter(result));
  CGAL::intersection(pgn_with_holes1, pgn_with_holes2,
                     std::back_inserter(result));
  CGAL::intersection(polygons.begin(), polygons.end(),
                     std::back_inserter(result));
  CGAL::intersection(polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result));
  CGAL::intersection(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result));

  CGAL::intersection(pgn1, pgn2, std::back_inserter(result), tr);
  CGAL::intersection(pgn1, pgn_with_holes2, std::back_inserter(result), tr);
  CGAL::intersection(pgn_with_holes1, pgn2, std::back_inserter(result), tr);
  CGAL::intersection(pgn_with_holes1, pgn_with_holes2,
                     std::back_inserter(result), tr);
  CGAL::intersection(polygons.begin(), polygons.end(),
                     std::back_inserter(result), tr);
  CGAL::intersection(polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result), tr);
  CGAL::intersection(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result), tr);

  Polygon_with_holes_2 res;
  CGAL::join(pgn1, pgn2, res);
  CGAL::join(pgn1, pgn_with_holes2, res);
  CGAL::join(pgn_with_holes1, pgn2, res);
  CGAL::join(pgn_with_holes1, pgn_with_holes2, res);
  CGAL::join(polygons.begin(), polygons.end(), std::back_inserter(result));
  CGAL::join(polygons_with_holes.begin(), polygons_with_holes.end(), std::back_inserter(result));
  CGAL::join(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(), std::back_inserter(result));

  CGAL::join(pgn1, pgn2, res, tr);
  CGAL::join(pgn1, pgn_with_holes2, res, tr);
  CGAL::join(pgn_with_holes1, pgn2, res, tr);
  CGAL::join(pgn_with_holes1, pgn_with_holes2, res, tr);
  CGAL::join(polygons.begin(), polygons.end(), std::back_inserter(result), tr);
  CGAL::join(polygons_with_holes.begin(), polygons_with_holes.end(),
             std::back_inserter(result), tr);
  CGAL::join(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(),
             std::back_inserter(result), tr);


  CGAL::difference(pgn1, pgn2, std::back_inserter(result));
  CGAL::difference(pgn1, pgn_with_holes2, std::back_inserter(result));
  CGAL::difference(pgn_with_holes1, pgn2, std::back_inserter(result));
  CGAL::difference(pgn_with_holes1, pgn_with_holes2,
                   std::back_inserter(result));

  CGAL::difference(pgn1, pgn2, std::back_inserter(result), tr);
  CGAL::difference(pgn1, pgn_with_holes2, std::back_inserter(result), tr);
  CGAL::difference(pgn_with_holes1, pgn2, std::back_inserter(result), tr);
  CGAL::difference(pgn_with_holes1, pgn_with_holes2,
                   std::back_inserter(result), tr);

  CGAL::symmetric_difference(pgn1, pgn2, std::back_inserter(result));
  CGAL::symmetric_difference(pgn1, pgn_with_holes2, std::back_inserter(result));
  CGAL::symmetric_difference(pgn_with_holes1, pgn2, std::back_inserter(result));
  CGAL::symmetric_difference(pgn_with_holes1, pgn_with_holes2,
                             std::back_inserter(result));
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             std::back_inserter(result));
  CGAL::symmetric_difference(polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result));
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result));

  CGAL::symmetric_difference(pgn1, pgn2, std::back_inserter(result), tr);
  CGAL::symmetric_difference(pgn1, pgn_with_holes2, std::back_inserter(result),
                             tr);
  CGAL::symmetric_difference(pgn_with_holes1, pgn2, std::back_inserter(result),
                             tr);
  CGAL::symmetric_difference(pgn_with_holes1, pgn_with_holes2,
                             std::back_inserter(result), tr);
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             std::back_inserter(result), tr);
  CGAL::symmetric_difference(polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result), tr);
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result), tr);

  Polygon_with_holes_2 res2;
  CGAL::complement(pgn1, res2);
  CGAL::complement(pgn_with_holes1, std::back_inserter(result));

  CGAL::complement(pgn1, res2, tr);
  CGAL::complement(pgn_with_holes1, std::back_inserter(result), tr);

  GPS gps2(pgn1);
  GPS gps3(pgn_with_holes1);

  GPS gps4;
  gps4.insert(pgn1);
  gps4.insert(pgn_with_holes1);
  gps4.insert(polygons.begin(), polygons.end());
  gps4.insert(polygons.begin(), polygons.end(),
              polygons_with_holes.begin(), polygons_with_holes.end());

  gps.complement(gps2);
  gps.intersection(gps2, gps3);
  gps.join(gps2, gps3);
  gps.difference(gps2, gps3);
  gps.symmetric_difference(gps2, gps3);

  Point_2 pt;
  gps.oriented_side(pt);
  gps.oriented_side(pgn1);
  gps.oriented_side(pgn_with_holes1);
  gps.oriented_side(pgn_with_holes2);
  gps.oriented_side(gps);
  gps.locate(pt, pgn_with_holes1);

  GPS new_gps(gps);
  GPS new_gps2 = gps;
}

void test_CGAL_Polygon_variants()
{
  typedef CGAL::Polygon_2<Kernel>               Polygon_2;
  typedef CGAL::Polygon_with_holes_2<Kernel>    Polygon_with_holes_2;
  typedef CGAL::Gps_default_traits<Polygon_2>::Traits Traits;

  Polygon_2 pgn1, pgn2;
  Polygon_with_holes_2  pgn_with_holes1, pgn_with_holes2;
  std::vector<Polygon_2>             polygons;
  std::vector<Polygon_with_holes_2>  polygons_with_holes;
  Polygon_with_holes_2 res;
  std::vector<Polygon_with_holes_2>  result;
  Traits tr;

  CGAL::do_intersect(pgn1, pgn2);
  CGAL::do_intersect(pgn1, pgn2, CGAL::Tag_true());
  CGAL::do_intersect(pgn1, pgn2, CGAL::Tag_false());
  CGAL::do_intersect(pgn1, pgn2, tr);

  CGAL::do_intersect(pgn1, pgn_with_holes2);
  CGAL::do_intersect(pgn1, pgn_with_holes2, CGAL::Tag_true());
  CGAL::do_intersect(pgn1, pgn_with_holes2, CGAL::Tag_false());
  CGAL::do_intersect(pgn1, pgn_with_holes2, tr);

  CGAL::do_intersect(pgn_with_holes1, pgn2);
  CGAL::do_intersect(pgn_with_holes1, pgn2, CGAL::Tag_true());
  CGAL::do_intersect(pgn_with_holes1, pgn2, CGAL::Tag_false());
  CGAL::do_intersect(pgn_with_holes1, pgn2, tr);

  CGAL::do_intersect(pgn_with_holes1, pgn_with_holes2);
  CGAL::do_intersect(pgn_with_holes1, pgn_with_holes2, CGAL::Tag_true());
  CGAL::do_intersect(pgn_with_holes1, pgn_with_holes2, CGAL::Tag_false());
  CGAL::do_intersect(pgn_with_holes1, pgn_with_holes2, tr);

  CGAL::do_intersect(polygons.begin(), polygons.end());
  CGAL::do_intersect(polygons.begin(), polygons.end(), CGAL::Tag_true());
  CGAL::do_intersect(polygons.begin(), polygons.end(), CGAL::Tag_false());
  CGAL::do_intersect(polygons.begin(), polygons.end(), tr);

  CGAL::do_intersect(polygons_with_holes.begin(), polygons_with_holes.end());
  CGAL::do_intersect(polygons_with_holes.begin(), polygons_with_holes.end(),
                     CGAL::Tag_true());
  CGAL::do_intersect(polygons_with_holes.begin(), polygons_with_holes.end(),
                     CGAL::Tag_false());
  CGAL::do_intersect(polygons_with_holes.begin(), polygons_with_holes.end(), tr);

  CGAL::do_intersect(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end());
  CGAL::do_intersect(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end(),
                     CGAL::Tag_true());
  CGAL::do_intersect(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end(),
                     CGAL::Tag_false());
  CGAL::do_intersect(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end(), tr);

  CGAL::intersection(pgn1, pgn2, std::back_inserter(result));
  CGAL::intersection(pgn1, pgn2, std::back_inserter(result), CGAL::Tag_true());
  CGAL::intersection(pgn1, pgn2, std::back_inserter(result), CGAL::Tag_false());
  CGAL::intersection(pgn1, pgn2, std::back_inserter(result), tr);

  CGAL::intersection(pgn1, pgn_with_holes2, std::back_inserter(result));
  CGAL::intersection(pgn1, pgn_with_holes2, std::back_inserter(result),
                     CGAL::Tag_true());
  CGAL::intersection(pgn1, pgn_with_holes2, std::back_inserter(result),
                     CGAL::Tag_false());
  CGAL::intersection(pgn1, pgn_with_holes2, std::back_inserter(result), tr);

  CGAL::intersection(pgn_with_holes1, pgn2, std::back_inserter(result));
  CGAL::intersection(pgn_with_holes1, pgn2, std::back_inserter(result),
                     CGAL::Tag_true());
  CGAL::intersection(pgn_with_holes1, pgn2, std::back_inserter(result),
                     CGAL::Tag_false());
  CGAL::intersection(pgn_with_holes1, pgn2, std::back_inserter(result), tr);

  CGAL::intersection(pgn_with_holes1, pgn_with_holes2,
                     std::back_inserter(result));
  CGAL::intersection(pgn_with_holes1, pgn_with_holes2,
                     std::back_inserter(result), CGAL::Tag_true());
  CGAL::intersection(pgn_with_holes1, pgn_with_holes2,
                     std::back_inserter(result), CGAL::Tag_false());
  CGAL::intersection(pgn_with_holes1, pgn_with_holes2,
                     std::back_inserter(result), tr);

  CGAL::intersection(polygons.begin(), polygons.end(),
                     std::back_inserter(result));
  CGAL::intersection(polygons.begin(), polygons.end(),
                     std::back_inserter(result), CGAL::Tag_true());
  CGAL::intersection(polygons.begin(), polygons.end(),
                     std::back_inserter(result), CGAL::Tag_false());
  CGAL::intersection(polygons.begin(), polygons.end(),
                     std::back_inserter(result), tr);

  CGAL::intersection(polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result));
  CGAL::intersection(polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result), CGAL::Tag_true());
  CGAL::intersection(polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result), CGAL::Tag_false());
  CGAL::intersection(polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result), tr);

  CGAL::intersection(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result));
  CGAL::intersection(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result), CGAL::Tag_true());
  CGAL::intersection(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result), CGAL::Tag_false());
  CGAL::intersection(polygons.begin(), polygons.end(),
                     polygons_with_holes.begin(), polygons_with_holes.end(),
                     std::back_inserter(result), tr);

  CGAL::join(pgn1, pgn2, res);
  CGAL::join(pgn1, pgn2, res, CGAL::Tag_true());
  CGAL::join(pgn1, pgn2, res, CGAL::Tag_false());
  CGAL::join(pgn1, pgn2, res, tr);

  CGAL::join(pgn1, pgn_with_holes2, res);
  CGAL::join(pgn1, pgn_with_holes2, res, CGAL::Tag_true());
  CGAL::join(pgn1, pgn_with_holes2, res, CGAL::Tag_false());
  CGAL::join(pgn1, pgn_with_holes2, res, tr);

  CGAL::join(pgn_with_holes1, pgn2, res);
  CGAL::join(pgn_with_holes1, pgn2, res, CGAL::Tag_true());
  CGAL::join(pgn_with_holes1, pgn2, res, CGAL::Tag_false());
  CGAL::join(pgn_with_holes1, pgn2, res, tr);

  CGAL::join(pgn_with_holes1, pgn_with_holes2, res);
  CGAL::join(pgn_with_holes1, pgn_with_holes2, res, CGAL::Tag_true());
  CGAL::join(pgn_with_holes1, pgn_with_holes2, res, CGAL::Tag_false());
  CGAL::join(pgn_with_holes1, pgn_with_holes2, res, tr);

  CGAL::join(polygons.begin(), polygons.end(), std::back_inserter(result));
  CGAL::join(polygons.begin(), polygons.end(), std::back_inserter(result),
             CGAL::Tag_true());
  CGAL::join(polygons.begin(), polygons.end(), std::back_inserter(result),
             CGAL::Tag_false());
  CGAL::join(polygons.begin(), polygons.end(), std::back_inserter(result), tr);

  CGAL::join(polygons_with_holes.begin(), polygons_with_holes.end(),
             std::back_inserter(result));
  CGAL::join(polygons_with_holes.begin(), polygons_with_holes.end(),
             std::back_inserter(result), CGAL::Tag_true());
  CGAL::join(polygons_with_holes.begin(), polygons_with_holes.end(),
             std::back_inserter(result), CGAL::Tag_false());
  CGAL::join(polygons_with_holes.begin(), polygons_with_holes.end(),
             std::back_inserter(result), tr);

  CGAL::join(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(),
             std::back_inserter(result));
  CGAL::join(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(),
             std::back_inserter(result), CGAL::Tag_true());
  CGAL::join(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(),
             std::back_inserter(result), CGAL::Tag_false());
  CGAL::join(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(),
             std::back_inserter(result), tr);

  CGAL::difference(pgn1, pgn2, std::back_inserter(result));
  CGAL::difference(pgn1, pgn2, std::back_inserter(result), CGAL::Tag_true());
  CGAL::difference(pgn1, pgn2, std::back_inserter(result), CGAL::Tag_false());
  CGAL::difference(pgn1, pgn2, std::back_inserter(result), tr);

  CGAL::difference(pgn1, pgn_with_holes2, std::back_inserter(result));
  CGAL::difference(pgn1, pgn_with_holes2, std::back_inserter(result),
                   CGAL::Tag_true());
  CGAL::difference(pgn1, pgn_with_holes2, std::back_inserter(result),
                   CGAL::Tag_false());
  CGAL::difference(pgn1, pgn_with_holes2, std::back_inserter(result), tr);

  CGAL::difference(pgn_with_holes1, pgn2, std::back_inserter(result));
  CGAL::difference(pgn_with_holes1, pgn2, std::back_inserter(result),
                   CGAL::Tag_true());
  CGAL::difference(pgn_with_holes1, pgn2, std::back_inserter(result),
                   CGAL::Tag_false());
  CGAL::difference(pgn_with_holes1, pgn2, std::back_inserter(result), tr);

  CGAL::difference(pgn_with_holes1, pgn_with_holes2, std::back_inserter(result));
  CGAL::difference(pgn_with_holes1, pgn_with_holes2, std::back_inserter(result),
                   CGAL::Tag_true());
  CGAL::difference(pgn_with_holes1, pgn_with_holes2, std::back_inserter(result),
                   CGAL::Tag_false());
  CGAL::difference(pgn_with_holes1, pgn_with_holes2, std::back_inserter(result),
                   tr);

  CGAL::symmetric_difference(pgn1, pgn2, std::back_inserter(result));
  CGAL::symmetric_difference(pgn1, pgn2, std::back_inserter(result),
                             CGAL::Tag_true());
  CGAL::symmetric_difference(pgn1, pgn2, std::back_inserter(result),
                             CGAL::Tag_false());
  CGAL::symmetric_difference(pgn1, pgn2, std::back_inserter(result), tr);

  CGAL::symmetric_difference(pgn1, pgn_with_holes2, std::back_inserter(result));
  CGAL::symmetric_difference(pgn1, pgn_with_holes2, std::back_inserter(result),
                             CGAL::Tag_true());
  CGAL::symmetric_difference(pgn1, pgn_with_holes2, std::back_inserter(result),
                             CGAL::Tag_false());
  CGAL::symmetric_difference(pgn1, pgn_with_holes2, std::back_inserter(result),
                             tr);

  CGAL::symmetric_difference(pgn_with_holes1, pgn2, std::back_inserter(result));
  CGAL::symmetric_difference(pgn_with_holes1, pgn2, std::back_inserter(result),
                             CGAL::Tag_true());
  CGAL::symmetric_difference(pgn_with_holes1, pgn2, std::back_inserter(result),
                             CGAL::Tag_false());
  CGAL::symmetric_difference(pgn_with_holes1, pgn2, std::back_inserter(result),
                             tr);

  CGAL::symmetric_difference(pgn_with_holes1, pgn_with_holes2,
                             std::back_inserter(result));
  CGAL::symmetric_difference(pgn_with_holes1, pgn_with_holes2,
                             std::back_inserter(result), CGAL::Tag_true());
  CGAL::symmetric_difference(pgn_with_holes1, pgn_with_holes2,
                             std::back_inserter(result), CGAL::Tag_false());
  CGAL::symmetric_difference(pgn_with_holes1, pgn_with_holes2,
                             std::back_inserter(result), tr);

  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             std::back_inserter(result));
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             std::back_inserter(result), CGAL::Tag_true());
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             std::back_inserter(result), CGAL::Tag_false());
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             std::back_inserter(result), tr);

  CGAL::symmetric_difference(polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result));
  CGAL::symmetric_difference(polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result), CGAL::Tag_true());
  CGAL::symmetric_difference(polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result), CGAL::Tag_false());
  CGAL::symmetric_difference(polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result), tr);

  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result));
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result), CGAL::Tag_true());
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result), CGAL::Tag_false());
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
                             polygons_with_holes.begin(),
                             polygons_with_holes.end(),
                             std::back_inserter(result), tr);

  Polygon_with_holes_2 res2;
  CGAL::complement(pgn1, res2);
  CGAL::complement(pgn1, res2, CGAL::Tag_true());
  CGAL::complement(pgn1, res2, CGAL::Tag_false());
  CGAL::complement(pgn1, res2, tr);

  CGAL::complement(pgn_with_holes1, std::back_inserter(result));
  CGAL::complement(pgn_with_holes1, std::back_inserter(result),
                   CGAL::Tag_true());
  CGAL::complement(pgn_with_holes1, std::back_inserter(result),
                   CGAL::Tag_false());
  CGAL::complement(pgn_with_holes1, std::back_inserter(result), tr);
}

int main()
{
  test<Gps>();
  test<Ps>();
  test<Gps_non_caching>();
  test_CGAL_Polygon_variants();

  return (0);
}
