// examples/Sweep_line/example5.C
// ------------------------------

#if defined(CGAL_CFG_NO_LONGNAME_PROBLEM) || defined(_MSC_VER)
#define Quotient                        Qt
#define Cartesian                       Cn
#define Arr_segment_exact_traits        AST
#define Pm_default_dcel                 PmDD
#define X_curve_plus_id                 XCPI
#define allocator                       All
#define Sweep_curves_to_subcurves_utils SCSU
#define Sweep_curves_base_2             SCB
#define Intersection_point_node         IPN
#endif

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/sweep_to_report_intersecting_curves_2.h>
#include <CGAL/IO/Pm_iostream.h>
#include <iostream>
#include <vector>

typedef CGAL::Quotient<CGAL::MP_Float>          NT;
typedef CGAL::Cartesian<NT>                     Kernel;
typedef CGAL::Arr_segment_exact_traits<Kernel>  Traits;

typedef Traits::Point_2                         Point_2;
typedef Traits::Curve_2                         Curve_2;
typedef Traits::X_curve_2                       X_curve_2;

typedef CGAL::Pm_default_dcel<Traits>           Dcel;   
typedef CGAL::Planar_map_2<Dcel,Traits>         PM;

typedef std::pair<Point_2, std::list<X_curve_2> > Pair_point_curves;

int main()
{
  PM                    pm;
  int                   num_segments;
  std::vector<Curve_2>  segments;
  
  std::cout << " * * * Demonstrating a trivial usage of the sweep line ";
  std::cout << "algorithm" << std::endl << std::endl;
  
  // Read input
  std::cin >> num_segments;
  
  while (num_segments--) 
  {
    Curve_2 cv;
    std::cin >> cv;
    segments.push_back(cv);
  }    

  // Use a sweep to find the points induced in the arrangement  
  std::vector<Pair_point_curves> intersecting_curves;
  Traits traits;
  CGAL::sweep_to_report_intersecting_curves_2(segments.begin(),
                                              segments.end(), 
                                              traits, 
                                              std::back_inserter(intersecting_curves));

  // Write output
  std::cout << " * * * Printing list of all points and their corresponding curves meeting at that point"<<std::endl;
  
  for (std::vector<Pair_point_curves>::iterator 
         p_iter = intersecting_curves.begin();
       p_iter != intersecting_curves.end(); ++p_iter)
    {
      std::cout<<"The curves"<< std::endl;
      for (std::list<X_curve_2>::iterator cv_iter = (p_iter->second).begin();
           cv_iter != (p_iter->second).end(); ++cv_iter)
        std::cout<< *cv_iter << std::endl;
      std::cout<<"Intersect in "<< p_iter->first << std::endl;
    }
  
  return 0;
}
