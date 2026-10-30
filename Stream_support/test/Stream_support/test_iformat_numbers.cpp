#include <CGAL/IO/io.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <type_traits>

int main()
{
  // 3D point coordinates
  {
    std::istringstream iss("1.5 2.25 3.75");
    double x, y, z;
    iss >> CGAL::IO::iformat(x) >> CGAL::IO::iformat(y) >> CGAL::IO::iformat(z);
    assert(iss && x == 1.5 && y == 2.25 && z == 3.75);
  }

  // float
  {
    std::istringstream iss("0.5 -2.0");
    float x, y;
    iss >> CGAL::IO::iformat(x) >> CGAL::IO::iformat(y);
    assert(iss && x == 0.5f && y == -2.0f);
  }

  // multi-line vertices with mixed whitespace
  {
    std::istringstream iss("\n 0 0 0\n\t1 0 0\n0 1 0\n0 0 1\n");
    double v[4][3];
    for(auto& p : v)
      iss >> CGAL::IO::iformat(p[0]) >> CGAL::IO::iformat(p[1]) >> CGAL::IO::iformat(p[2]);
    assert(iss);
    assert(v[1][0] == 1 && v[2][1] == 1 && v[3][2] == 1 && v[3][0] == 0);
  }

  // integers, including limits
  {
    std::istringstream iss("0 1 2 2147483647 -2147483648");
    int i[5];
    for(int& k : i)
      iss >> CGAL::IO::iformat(k);
    assert(iss);
    assert(i[0] == 0 && i[1] == 1 && i[2] == 2 && i[3] == 2147483647 && i[4] == -2147483647 - 1);
  }

  // the delimiter after a number must remain in the stream
  {
    std::istringstream iss("42, 100");
    int x, y;
    char comma;
    iss >> CGAL::IO::iformat(x) >> comma >> CGAL::IO::iformat(y);
    assert(iss && x == 42 && comma == ',' && y == 100);
  }

  // reading stops at the first character that cannot be part of a number
  {
    std::istringstream iss("42,100;-1.5e-3)2.5");
    int x, y;
    double d, e;
    char c1, c2, c3;
    iss >> CGAL::IO::iformat(x) >> c1 >> CGAL::IO::iformat(y) >> c2
        >> CGAL::IO::iformat(d) >> c3 >> CGAL::IO::iformat(e);
    assert(iss && x == 42 && c1 == ',' && y == 100 && c2 == ';');
    assert(d == -1.5e-3 && c3 == ')' && e == 2.5);
  }

  // other integral types
  {
    std::istringstream iss("-5 7 4000000000 -9000000000 65535");
    short s; unsigned u; unsigned int ui; long long ll; unsigned short us;
    iss >> CGAL::IO::iformat(s) >> CGAL::IO::iformat(u) >> CGAL::IO::iformat(ui)
        >> CGAL::IO::iformat(ll) >> CGAL::IO::iformat(us);
    assert(iss && s == -5 && u == 7 && ui == 4000000000u && ll == -9000000000LL && us == 65535);
  }

  // out-of-range values and a lone sign fail
  {
    std::istringstream iss("3000000000");
    int x = 1;
    iss >> CGAL::IO::iformat(x);
    assert(iss.fail());
  }
  {
    std::istringstream iss("70000");
    unsigned short x = 1;
    iss >> CGAL::IO::iformat(x);
    assert(iss.fail());
  }
  {
    std::istringstream iss("1e999");
    double x = 1;
    iss >> CGAL::IO::iformat(x);
    assert(iss.fail());
  }
  {
    std::istringstream iss("- 5");
    int x = 1;
    iss >> CGAL::IO::iformat(x);
    assert(iss.fail());
  }

  // long double is not handled by from_chars and uses the generic operator>>
  {
    std::istringstream iss("1.5 2");
    long double x, y;
    iss >> CGAL::IO::iformat(x) >> CGAL::IO::iformat(y);
    assert(iss && x == 1.5L && y == 2);
    static_assert(std::is_base_of_v<CGAL::IO_rep_is_not_specialized, CGAL::Input_rep<long double> >);
  }

  // char and bool are still read as characters, not numbers
  {
    std::istringstream iss("a 1");
    char c; bool b;
    iss >> CGAL::IO::iformat(c) >> CGAL::IO::iformat(b);
    assert(iss && c == 'a' && b);
    static_assert(std::is_base_of_v<CGAL::IO_rep_is_not_specialized, CGAL::Input_rep<char> >);
    static_assert(std::is_base_of_v<CGAL::IO_rep_is_not_specialized, CGAL::Input_rep<bool> >);
#ifndef CGAL_DISABLE_FROM_CHARS
    static_assert(std::is_base_of_v<CGAL::IO_rep_is_specialized, CGAL::Input_rep<unsigned long> >);
#endif
  }

  // last number without trailing whitespace
  {
    std::istringstream iss("999");
    int x;
    iss >> CGAL::IO::iformat(x);
    assert(iss && x == 999);
  }

  // empty/blank stream fails and leaves the value untouched
  {
    std::istringstream iss("   ");
    int x = 999;
    iss >> CGAL::IO::iformat(x);
    assert(iss.fail() && x == 999);
  }

  // invalid input
  {
    std::istringstream iss("123 invalid 456");
    int x;
    double y = 7;
    iss >> CGAL::IO::iformat(x);
    assert(iss && x == 123);
    iss >> CGAL::IO::iformat(y);
    assert(iss.fail());
  }

  std::cout << "done" << std::endl;
  return 0;
}
