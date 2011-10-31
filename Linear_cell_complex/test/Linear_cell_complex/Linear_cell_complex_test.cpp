#include "Linear_cell_complex_2_test.h"
#include "Linear_cell_complex_3_test.h"
#include "Linear_cell_complex_4_test.h"

int main()
{
  typedef CGAL::Linear_cell_complex<2> LCC1;
  if ( !test_LCC_2<LCC1>() )
  {
    std::cout<<"ERROR during Test_LCC_2<LCC1>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
