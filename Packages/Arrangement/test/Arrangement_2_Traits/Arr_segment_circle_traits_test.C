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

#include <CGAL/leda_real.h>
#include "include/Segment_circle_traits_test.h"
#include <CGAL/Arr_segment_circle_traits.h>

typedef leda_real                                  NT;
typedef CGAL::Arr_segment_circle_traits<NT>        Traits;

int main (int argc, char** argv)
{
  Segment_circle_traits_test<Traits, NT>  test_obj (argc, argv);

  if (test_obj.start())
    return (0); // SUCCESS
  else 
    return (1); // FAILURE  
}

#endif // ! defined(CGAL_USE_LEDA) ...
