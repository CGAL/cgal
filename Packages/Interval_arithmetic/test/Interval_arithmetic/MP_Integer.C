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

int main()
{
  MPI a(0);
  for (int i=0; i<10000; i++)
    a = a + MPI(i);

  std::cout << a << std::endl;
  CGAL_assertion ( a == MPI(5000*9999));
  
  MPI b = factoriel(10);
  std::cout << b << std::endl;
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
