// provide 3d kernel traits ...
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3 
 
#include <CGAL/basic.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include "CGAL/_test_new_2.h"

typedef CGAL::leda_rat_kernel_traits   Kernel;

int main()
{
  test_new_2( Kernel() );
  return 0;
}

