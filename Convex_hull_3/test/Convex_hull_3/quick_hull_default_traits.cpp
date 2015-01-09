#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>

#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>


typedef CGAL::Exact_rational Exact_rational;


typedef CGAL::Exact_predicates_exact_constructions_kernel   EPEC;
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;
typedef CGAL::Simple_cartesian<double>                      SCD;
typedef CGAL::Simple_cartesian<Exact_rational>            SCR;

using namespace CGAL::internal::Convex_hull_3;

int main()
{
  CGAL_static_assertion( (boost::is_same<EPEC,Default_traits_for_Chull_3<EPEC::Point_3>::type>::value) );
  CGAL_static_assertion( (boost::is_same<SCD,Default_traits_for_Chull_3<SCD::Point_3>::type>::value) );
  CGAL_static_assertion( (boost::is_same<SCR,Default_traits_for_Chull_3<SCR::Point_3>::type>::value) );
  CGAL_static_assertion( (boost::is_same<EPEC,Default_traits_for_Chull_3<EPEC::Point_3>::type>::value) );
  CGAL_static_assertion( (boost::is_same<CGAL::Convex_hull_traits_3<EPIC, CGAL::Tag_true>,Default_traits_for_Chull_3<EPIC::Point_3>::type>::value) );
  CGAL_static_assertion( (boost::is_same<Is_on_positive_side_of_plane_3<CGAL::Convex_hull_traits_3<EPIC, CGAL::Tag_true> >::Protector,CGAL::Protect_FPU_rounding<true> >::value) );
  return 0;
}
