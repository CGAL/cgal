#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_rational.h>
typedef leda_rational                   Precise_rational;
#elif defined CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#  include <CGAL/Quotient.h>
typedef CGAL::Quotient<CGAL::Gmpz>      Precise_rational;
#else
#  include <CGAL/MP_Float.h>
#  include <CGAL/Quotient.h>
typedef CGAL::Quotient<CGAL::MP_Float>  Precise_rational;
#endif

typedef CGAL::Exact_predicates_exact_constructions_kernel   EPEC;
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;
typedef CGAL::Simple_cartesian<double>                      SCD;
typedef CGAL::Simple_cartesian<Precise_rational>            SCR;

using namespace CGAL::internal::Convex_hull_3;

int main()
{
  BOOST_STATIC_ASSERT( (boost::is_same<EPEC,Default_traits_for_Chull_3<EPEC::Point_3>::type>::value) );
  BOOST_STATIC_ASSERT( (boost::is_same<SCD,Default_traits_for_Chull_3<SCD::Point_3>::type>::value) );
  BOOST_STATIC_ASSERT( (boost::is_same<SCR,Default_traits_for_Chull_3<SCR::Point_3>::type>::value) );
  BOOST_STATIC_ASSERT( (boost::is_same<EPEC,Default_traits_for_Chull_3<EPEC::Point_3>::type>::value) );
  BOOST_STATIC_ASSERT( (boost::is_same<CGAL::Convex_hull_traits_3<EPIC>,Default_traits_for_Chull_3<EPIC::Point_3>::type>::value) );
  BOOST_STATIC_ASSERT( (boost::is_same<Is_on_positive_side_of_plane_3<CGAL::Convex_hull_traits_3<EPIC> >::Protector,CGAL::Protect_FPU_rounding<true> >::value) );
  return 0;
}
