#include <CGAL/Exact_rational.h>
#include <iostream>

typedef CGAL::Exact_rational Rational;

int main()
{
  std::cout << typeid(Rational).name() << std::endl;

  Rational rat(-0.375);
  double d = -0.625;

  std::cout << "rat : " << rat << std::endl;

  Rational sub = rat - d;
  std::cout << "sub : " << sub << std::endl;
  assert(rat != sub);

  return 0;
}
