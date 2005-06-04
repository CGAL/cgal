#ifndef SVD_MAIN_HIERARCHY_H
#define SVD_MAIN_HIERARCHY_H

#include <CGAL/basic.h>
#include <iostream>
#include <cassert>

#include "IO/io_aux.h"

bool assert_no_warning(bool b)
{
  assert(b);
  return b;
}

int main(int argc, char* argv[])
{
  {
    print_separator();

    SVD2 svd;

    bool types_ok = CGAL::test_svd_hierarchy(std::cin, svd, "svdh.tmp");

    assert_no_warning( types_ok );

    print_separator();

    std::cout << std::endl;
  }

  {
    print_separator();

    SVD2_wi svd;

    bool types_ok = CGAL::test_svd_hierarchy(std::cin, svd, "svdh_wi.tmp");

    assert_no_warning( types_ok );

    print_separator();

    std::cout << std::endl;
  }

  return 0;
}

#endif // SVD_MAIN_HIERARCHY_H
