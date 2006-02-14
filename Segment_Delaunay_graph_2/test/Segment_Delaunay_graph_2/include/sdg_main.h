#ifndef SDG_MAIN_H
#define SDG_MAIN_H

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

    bool types_ok = CGAL::test_sdg(std::cin, sdg, "sdg.tmp");

    assert_no_warning( types_ok );

    print_separator();

    std::cout << std::endl;
  }

  {
    print_separator();

    SDG2_wi sdg;

    bool types_ok = CGAL::test_sdg(std::cin, sdg, "sdg_wi.tmp");

    assert_no_warning( types_ok );

    print_separator();

    std::cout << std::endl;
  }

  return 0;
}

#endif // SDG_MAIN_H
