#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_segment_cached_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>

#include "include/Polyline_traits_test.h"

typedef CGAL::Quotient<CGAL::MP_Float>             NT;
typedef CGAL::Cartesian<NT>                        Kernel;
typedef CGAL::Arr_segment_cached_traits_2<Kernel>  Seg_traits;
typedef CGAL::Arr_polyline_traits_2<Seg_traits>    Traits;

int main( int argc, char** argv )
{
  Polyline_traits_test< Traits, NT >  test_obj( argc, argv );

  if (test_obj.start())
    return (0); // SUCCESS
  else 
    return (1); // FAILURE  
}
