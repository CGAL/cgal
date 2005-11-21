#ifndef SDG_MAIN_HIERARCHY_H
#define SDG_MAIN_HIERARCHY_H

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

    SDG2 sdg;

    bool types_ok = CGAL::test_sdg_hierarchy(std::cin, sdg, "sdgh.tmp");

    assert_no_warning( types_ok );

    print_separator();

    std::cout << std::endl;
  }

  {
    print_separator();

    SDG2_wi sdg;

    bool types_ok = CGAL::test_sdg_hierarchy(std::cin, sdg, "sdgh_wi.tmp");

    assert_no_warning( types_ok );

    print_separator();

    std::cout << std::endl;
  }

  return 0;
}

#endif // SDG_MAIN_HIERARCHY_H
