#include <CGAL/basic.h>

#define CGAL_NO_DEPRECATED_CODE
#define CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS

#include <CGAL/Kernel_archetype.h>

// needed in kernel testsuite ...
CGAL::Test_vector_2 operator-(const CGAL::Test_vector_2& v)
{ return v; }  

#include "CGAL/_test_new_2.h"

typedef CGAL::Kernel_archetype   Kernel;

int main()
{
  test_new_2( Kernel() );
  return 0;
}

