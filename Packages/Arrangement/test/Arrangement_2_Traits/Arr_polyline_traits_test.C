#include <CGAL/Cartesian.h>

#include <CGAL/Arr_polyline_traits.h>
#include "include/Polyline_traits_test.h"
#include <CGAL/Quotient.h>

typedef CGAL::Quotient<int>            NT;
typedef CGAL::Cartesian<NT>            Rep;
typedef CGAL::Arr_polyline_traits<Rep> Traits;

int main( int argc, char** argv ){
  Polyline_traits_test< Traits, NT >  test_obj( argc, argv );

  if (test_obj.start())
    return (0); // SUCCESS
  else 
    return (1); // FAILURE  
}
