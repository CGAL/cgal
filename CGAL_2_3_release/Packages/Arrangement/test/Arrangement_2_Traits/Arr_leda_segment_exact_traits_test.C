#include <CGAL/basic.h>
#include <iostream>

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA)
int main(int argc, char* argv[])
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else

#include <CGAL/Arr_leda_segment_exact_traits.h>
#include <CGAL/leda_rational.h>

#include "include/Segment_traits_test.h"

typedef leda_rational                                  NT;
typedef CGAL::Arr_leda_segment_exact_traits            Traits;

int main( int argc, char** argv ){
  Segment_traits_test< Traits, leda_rational >  test_obj( argc, argv );

  if (test_obj.start())
    return (0); // SUCCESS
  else 
    return (1); // FAILURE  
}

#endif // ! defined(CGAL_USE_LEDA) ...
