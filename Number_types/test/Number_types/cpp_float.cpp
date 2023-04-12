
#include <CGAL/cpp_float.h>

int main()
{

  CGAL::cpp_float m(0);
  std::cout << m << std::endl;

  CGAL::cpp_float m0(23.0);

  CGAL::cpp_float m1(5.0);

  CGAL::cpp_float m2(5.125);

  CGAL::cpp_float m3(2.5);

  CGAL::cpp_float m4(0.625);


  CGAL::cpp_float m5(0.5);

  CGAL::is_positive(m5);

  assert(m4 > m5);

  assert(-m4 < -m5);


  assert(-m4 == -m4);


  assert(-m4 != -m5);

  assert(-m5 != -m4);

  CGAL::cpp_float m15 = m1 + m5;
  m15 = m15 - m5;
  assert(m15 == m1);
  assert(! (m15 < m1));
  assert(! (m1 < m15));

  m15 = m15 - m15;
  std::cout << m15 << std::endl;
  assert(m15.is_zero());

  m1 *= m5;

  m0 += m4;
  std::cout << m0 << std::endl;

  return 0;
}
