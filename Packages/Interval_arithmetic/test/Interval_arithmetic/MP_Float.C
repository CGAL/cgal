// Test program for the MP_Integer class.
// Sylvain Pion.
//
// We should be able to more deeply test it using an NT checker.

#include <CGAL/MP_Integer.h>
#include <iostream>

typedef CGAL::MP_Integer MPI;

void test_equality(int i)
{
  assert(MPI(i) == MPI(i));
  assert(MPI(i) <= MPI(i));
  assert(MPI(i) >= MPI(i));
  assert(! (MPI(i) < MPI(i)));
  assert(! (MPI(i) > MPI(i)));
}

void bench(int loops)
{
  MPI j (13579);
  MPI k (9753);
  MPI l (-149);
  for(; loops>0; loops--)
  {
    MPI m (loops);
    MPI n = (j*k-l*m)*(j*m-k*l)*(j*j);
    (void) (n>MPI(0));
  }
}

MPI factoriel (short i)
{
  MPI r(1);
  for (short j=1; j<=i; j++)
    r = r * MPI(j);
  return r;
}

void print_test()
{
  for (int i=-2; i<3; i++)
    std::cout << MPI(i) << "   ";
  std::cout << endl;
  for (int i=-2; i<3; i++)
    std::cout << MPI(i+65536) << "      ";
  std::cout << endl;
}

int main(int argc, char **argv)
{
#if 0
  std::vector<int> iii;
  iii.push_back(1);
  iii.push_back(2);
  iii.push_back(3);
  iii.erase(iii.begin(), iii.begin()+1);
  std::cout << *iii.begin();
#endif

  int loops = argc > 1 ? atoi(argv[1]) : 0;
  bench(loops);

  std::cout.precision(20);

  for (; loops >= 0; loops--) {
    double d = drand48();
    int exp = int((drand48()-.5)*1024);
    d = std::ldexp(d, exp);
    // d = 1.9696110926449043849e+124;
    // d = 7.3074228557478900057e+47;
    // std::cout << d << std::endl;
    // std::cout << MPI(d) << std::endl;
    // std::cout << CGAL::to_double(MPI(d)) << std::endl;
    assert(CGAL::to_double(MPI(d)) == d);
    // MPI z(d);
  }

  MPI a(0);
  int j=0;
  for (int i=0; i<10000; i++)
  {
#if 0
    std::cout << i << std::endl;
    std::cout << a << std::endl;
    std::cout << "exp = " << a.exp << "    v[" << a.v.size()<< "] = ";
    if (!a.is_zero())
      for (int j=0; j<a.v.size(); j++)
        std::cout << " " << a.v[j];
    std::cout << std::endl;
#endif
    a = a + MPI(i);
    j += i;
    CGAL_assertion(a == MPI(j));
  }

  std::cout << a << std::endl;
  CGAL_assertion ( a == MPI(5000*9999));

  MPI bb = factoriel(100);
  std::cout << "100! = " << bb << std::endl;
  
  MPI b = factoriel(10);
  std::cout << "10! = " << b << " =? 3628800 " << " =? " << CGAL::to_double(b);
  std::cout << std::endl;
  CGAL_assertion ( b == MPI(3628800));

  test_equality(0);
  test_equality(10);
  test_equality(-10);

  MPI c = -10*65536+30; // Note : limb==short is hardcoded here
  MPI d = c * 65536;

  CGAL_assertion(c < (c+1));
  CGAL_assertion(c > (c-1));
  CGAL_assertion(d < (d+1));
  CGAL_assertion(d > (d-1));

  MPI e =  7*65536 + 10;
  MPI f = 10*65536 + 7;
  CGAL_assertion(e < f);

  print_test();

  return 0;
}
