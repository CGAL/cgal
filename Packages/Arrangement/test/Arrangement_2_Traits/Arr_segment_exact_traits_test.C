#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Quotient.h>

#include "include/Segment_traits_test.h"

typedef CGAL::Quotient<int>                            NT;
typedef CGAL::Cartesian<NT>                            R;
typedef CGAL::Arr_segment_traits_2<R>                  Traits;

int main(int argc, char * argv[])
{
  Segment_traits_test< Traits, NT >  test_obj( argc, argv );

  if (test_obj.start())
    return (0); // SUCCESS    
  else 
    return (1); // FAILURE  
}
