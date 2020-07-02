#include <CGAL/config.h>
#include <CGAL/use.h>
#include <iostream>

#if (!CGAL_USE_CORE) && (!CGAL_USE_LEDA) && (!(CGAL_USE_GMP && CGAL_USE_MPFI))
int main ()
{
  std::cout << "Sorry, this example needs CORE, LEDA, or GMP+MPFI ..."
            << std::endl;
  return 0;
}
#else

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#if CGAL_USE_GMP && CGAL_USE_MPFI
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz Integer;
#elif CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
typedef CORE::BigInt Integer;
#else
#include <CGAL/leda_integer.h>
typedef LEDA::integer Integer;
#endif

typedef CGAL::Arr_algebraic_segment_traits_2<Integer>   Arr_traits_2;
typedef CGAL::Arrangement_2<Arr_traits_2>               Arrangement_2;
typedef Arr_traits_2::Curve_2                           Curve_2;
typedef Arr_traits_2::Polynomial_2                      Polynomial_2;
typedef Arr_traits_2::Algebraic_real_1                  Algebraic_real_1;
typedef Arr_traits_2::X_monotone_curve_2                X_monotone_curve_2;
typedef Arr_traits_2::Point_2                           Point_2;

typedef boost::variant<Point_2, X_monotone_curve_2>     Make_x_monotone_result;

int main()
{
  Arr_traits_2 arr_traits;
  auto construct_curve = arr_traits.construct_curve_2_object();
  auto construct_x_monotone_segment =
    arr_traits.construct_x_monotone_segment_2_object();
  auto construct_point = arr_traits.construct_point_2_object();
  auto make_x_monotone = arr_traits.make_x_monotone_2_object();

  Arrangement_2 arr(&arr_traits);

  std::vector<X_monotone_curve_2> segs;

  Polynomial_2 x = CGAL::shift(Polynomial_2(1),1,0);
  Polynomial_2 y = CGAL::shift(Polynomial_2(1),1,1);

  // Construct x^4+y^3-1
  Curve_2 cv0 = construct_curve(CGAL::ipower(x,4)+CGAL::ipower(y,3)-1);
  // Construct all x-monotone segments using the Make_x_mononotone functor
  std::vector<Make_x_monotone_result> pre_segs;
  make_x_monotone(cv0, std::back_inserter(pre_segs));
  // Cast all CGAL::Objects into X_monotone_segment_2
  // (the vector might also contain Point_2 objects for isolated points,
  // but not for this instance
  for(size_t i = 0; i < pre_segs.size(); ++i) {
    auto* curr_p = boost::get<X_monotone_curve_2>(&pre_segs[i]);;
    CGAL_assertion(curr_p);
    segs.push_back(*curr_p);
  }
  // Construct an ellipse with equation 2*x^2+5*y^2-7=0
  Curve_2 cv1 = construct_curve(2*CGAL::ipower(x,2)+5*CGAL::ipower(y,2)-7);

  // Construct point on the upper arc (counting of arc numbers starts with 0!
  Point_2 p11 = construct_point(Algebraic_real_1(0),cv1,1);

  construct_x_monotone_segment(cv1,p11,Arr_traits_2::POINT_IN_INTERIOR,
                               std::back_inserter(segs));

  // Construct a vertical cusp x^2-y^3=0
  Curve_2 cv2 = construct_curve(CGAL::ipower(x,2)-CGAL::ipower(y,3));

  // Construct a segment containing the cusp point.
  // This adds to X_monotone_curve_2 objects to the vector,
  // because the cusp is a critical point
  Point_2 p21 = construct_point(Algebraic_real_1(-2),cv2,0);
  Point_2 p22 = construct_point(Algebraic_real_1(2),cv2,0);
  construct_x_monotone_segment(cv2,p21,p22,std::back_inserter(segs));

  // Construct an unbounded curve, starting at x=3
  Point_2 p23 = construct_point(Algebraic_real_1(3),cv2,0);
  construct_x_monotone_segment(cv2,p23,Arr_traits_2::MIN_ENDPOINT,
                                 std::back_inserter(segs));

  // Construct another conic: y^2-x^2+1
  Curve_2 cv3 = construct_curve(CGAL::ipower(y,2)-CGAL::ipower(x,2)+1);

  Point_2 p31 = construct_point(Algebraic_real_1(2),cv3,1);
  construct_x_monotone_segment(cv3,p31,Arr_traits_2::MAX_ENDPOINT,
                               std::back_inserter(segs));

  // Construct a vertical segment
  Point_2 v1 = construct_point(0,0);
  Point_2 v2 = construct_point(Algebraic_real_1(0),cv1,1);
  construct_x_monotone_segment(v1,v2,std::back_inserter(segs));

  CGAL::insert(arr,segs.begin(),segs.end());

  // Add some isolated points (must be wrapped into CGAL::Object)
  std::vector<CGAL::Object> isolated_points;
  isolated_points.push_back
    (CGAL::make_object(construct_point(Algebraic_real_1(2),cv3,0)));
  isolated_points.push_back
    (CGAL::make_object(construct_point(Integer(1),Integer(5))));
  isolated_points.push_back
    (CGAL::make_object(construct_point(Algebraic_real_1(-1),
                                       Algebraic_real_1(5))));

  CGAL::insert(arr,isolated_points.begin(), isolated_points.end());

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  return 0;
}

#endif
