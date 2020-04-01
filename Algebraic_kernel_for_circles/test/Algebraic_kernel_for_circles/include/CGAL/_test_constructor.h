#include <CGAL/Random.h>
#include <cassert>

template <class AK>
void _test_constuctor(AK ak)
{
  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127;
  typename AK::Construct_polynomial_for_circles_2_2 theConstruct_2_2 =
    ak.construct_polynomial_for_circles_2_2_object();
  typename AK::Construct_polynomial_1_2 theConstruct_1_2 =
    ak.construct_polynomial_1_2_object();

  for(int i = 0; i < 20 ; i++){
    int x = theRandom.get_int(random_min,random_max);
    int y = theRandom.get_int(random_min,random_max);
    int r_sq = theRandom.get_int(random_min,random_max);
    int a = theRandom.get_int(random_min,random_max);
    int b = theRandom.get_int(random_min,random_max);
    int c = theRandom.get_int(random_min,random_max);

    typename AK::Polynomial_for_circles_2_2 p_2_2 = theConstruct_2_2(x, y, r_sq);
    typename AK::Polynomial_1_2 p_1_2 = theConstruct_1_2(a, b, c);
    assert(p_2_2.a() == x);
    assert(p_2_2.b() == y);
    assert(p_2_2.r_sq() == r_sq);
    assert(p_1_2.a() == a);
    assert(p_1_2.b() == b);
    assert(p_1_2.c() == c);
  }

}
