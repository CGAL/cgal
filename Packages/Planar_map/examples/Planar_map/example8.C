/*! \file
 * A construction of a house-like-shaped Planar_map.
 * In this example a house-like-shaped Planar_map that uses the naive point
 * location strategy is constructed and printed out to the standard output.
 * By default, the planar map is templated with the straight traits, which is
 * templated with the cartesian kernel over doubles in turn. 
 * It is possible to toggle between naive and default point location strategy,
 * between cartesian and homogeneous kernel, and between double and leda
 * rational number type.
 */

#define CGAL_NO_PM_DEFAULT_POINT_LOCATION
#define NAIVE_POINT_LOCATION
#define CARTESIAN

#if defined(CGAL_USE_LEDA)
#include <CGAL/leda_rational.h>
#endif

#ifdef CARTESIAN
#include <CGAL/Quotient.h>
#include <CGAL/Cartesian.h>
#else
#include <CGAL/Homogeneous.h>
#endif

#include <CGAL/Pm_straight_traits.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/IO/Straight_2_stream.h>
#include <CGAL/IO/Pm_file_writer.h>
#include <CGAL/IO/write_pm.h>

#ifdef NAIVE_POINT_LOCATION
#include <CGAL/Pm_naive_point_location.h>
#endif

#include <iostream>
#include <iterator>
#include <algorithm>

#if defined(CGAL_USE_LEDA) && defined(USE_LEDA)
#ifdef CARTESIAN
typedef CGAL::Cartesian<leda_rational>          Kernel;
#else
typedef CGAL::Homogeneous<leda_rational>        Kernel;
#endif
#else //defined(CGAL_USE_LEDA) && defined(USE_LEDA)
#ifdef CARTESIAN
typedef CGAL::Quotient<float>                   Number_type;
typedef CGAL::Cartesian<Number_type>            Kernel;
#else
typedef CGAL::Homogeneous<double>               Kernel;
#endif
#endif

typedef CGAL::Pm_straight_traits<Kernel>        Traits;
typedef Traits::Bounding_box                    Bounding_box;
typedef Traits::Point_2                         Point_2;
typedef Traits::X_curve_2                       X_curve_2;
typedef Traits::X_bounded_curve                 Segment;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;
typedef Planar_map::Halfedge_handle             Halfedge_handle;
typedef CGAL::Pm_file_writer<Planar_map>        Pm_writer;

CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(X_curve_2);

int main()
{
  Bounding_box bbox(Point_2(-10, -10), Point_2(10, 10));
  Traits traits(bbox);

  // Create an instance of a planar map:
#ifndef NAIVE_POINT_LOCATION   
  Planar_map pm(traits, NULL, NULL);
#else
  CGAL::Pm_naive_point_location< Planar_map > naive_pl;
  Planar_map pm(traits, &naive_pl, NULL);
#endif
  Pm_writer verbose_writer(std::cout, pm, true);
  
  X_curve_2 cv[6];
  int i;
  
  CGAL::set_ascii_mode(std::cout);
  Point_2 a1(1, 1), a2(1, 0), a3(0, 0), a4(0, 1), a5(1, 4, 2) ;
  
  /*
        a5 
        /\  
       /  \
   a4  ----  a1
      |    | 
      |    | 
      |    |
   a3  ----  a2
 */        

  // Create the curves:
  cv[0] = X_curve_2(Segment(a1, a2));
  cv[1] = X_curve_2(Segment(a2, a3));
  cv[2] = X_curve_2(Segment(a3, a4));
  cv[3] = X_curve_2(Segment(a4, a5));
  cv[4] = X_curve_2(Segment(a5, a1));
  cv[5] = X_curve_2(Segment(a1, a4));
  
  std::cout << "The curves of the map :" << std::endl; 
  std::copy(&cv[0], &cv[6], std::ostream_iterator<X_curve_2>(std::cout, "\n"));
  std::cout << std::endl;
  
  // insert the five curves to the map
  std::cout << "Inserting the curves into the map ... ";
  Halfedge_handle e[6];
  e[0] = pm.insert(cv[0]);
  for (i = 1; i < 4; i++)
    e[i] = pm.insert_from_vertex(cv[i], e[i-1]->target());

  e[4] = pm.insert_at_vertices(cv[4], e[0]->source(), e[3]->target());
  e[5] = pm.insert_at_vertices(cv[5], e[0]->source(), e[2]->target());
  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;

  /* 
         O            
     e3 /\  e4
       /  \
      O----O
      | e5 | 
   e2 |    | e0
      |    |
      O----O
        e1 
 */

  CGAL::write_pm(pm, verbose_writer, std::cout);

  return 0;  
}
