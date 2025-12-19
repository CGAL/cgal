#ifdef CGAL_DO_NOT_USE_BOOST_MP

#include <iostream>

int main()
{
  std::cout << "The class CGAL::cpp_float is not tested on this platform" << std::endl;
  return 0;
}

#else

#include <CGAL/cpp_float.h>

template<class NT>
void test1(){
  NT z;
  NT a=3;
  NT b=4.5;
  NT c=2*(a+b)+-a*5;
  assert(CGAL::sign(c)==0);

  NT e=.0003;
  NT f=1e-90;
  assert(CGAL::to_double(b) == 4.5);
  std::pair<double,double> p=CGAL::to_interval(b);
  assert(p.first<=4.5 && p.second >= 4.5);
  assert(a<b && CGAL::compare(b,a)>0);
  assert(z<f && CGAL::compare(z,f)<0 && CGAL::compare(f,z)>0);
  assert(z==z && CGAL::compare(z,z)==0);
  assert(CGAL::square(b)*4==81);
  assert(CGAL::is_zero(c));
  assert(!CGAL::is_zero(a));
  assert(!CGAL::is_one(a));
  assert(CGAL::is_one(a-2));
  assert(e-e==0);
}

int main()
{
  test1<CGAL::cpp_float>();

  double d = -0;
  CGAL::cpp_float zero(d);
  assert(! CGAL::is_positive(zero));
  assert(! CGAL::is_negative(zero));
  assert(CGAL::is_zero(zero));

  CGAL::cpp_float m(0);
  std::cout << m << std::endl;

  CGAL::cpp_float m0(23.0);

  CGAL::cpp_float m1(5.0);

  CGAL::cpp_float m2(5.125);

  CGAL::cpp_float m3(2.5);

  CGAL::cpp_float m4(0.625);


  CGAL::cpp_float m5(0.5);

  CGAL::cpp_float m6(-0.625);

  CGAL::is_positive(m5);
  CGAL::is_negative(m6);

  assert(m < m5);
  assert(m6 < m);
  assert(! (m5 < m));
  assert(! (m < m6));
  assert(! (m6 < m6));



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

  CGAL::cpp_float one(1);
  assert(one.is_one());
  one += m4;
  one -= m4;
  std::cout << one << std::endl;
  assert(one.is_one());


  return 0;
}

#endif
