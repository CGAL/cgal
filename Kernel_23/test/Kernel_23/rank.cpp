#include <CGAL/rank.h>
#include <cassert>
#include <iostream>

typedef int FT;


void test_rank_33(FT a0, FT b0, FT c0,
                  FT a1, FT b1, FT c1,
                  FT a2, FT b2, FT c2,
                  int expected)
{
  std::cout << "testing:\n"
            << "\t " << a0 << "\t" << b0 << "\t" << c0 << "\n"
            << "\t " << a1 << "\t" << b1 << "\t" << c1 << "\n"
            << "\t " << a2 << "\t" << b2 << "\t" << c2 << "\n";

  assert(CGAL::rank_33<FT>(a0,b0,c0,a1,b1,c1,a2,b2,c2)==expected);
  assert(CGAL::rank_33<FT>(a0,b0,c0,a2,b2,c2,a1,b1,c1)==expected);
  assert(CGAL::rank_33<FT>(a1,b1,c1,a0,b0,c0,a2,b2,c2)==expected);
  assert(CGAL::rank_33<FT>(a1,b1,c1,a2,b2,c2,a0,b0,c0)==expected);
  assert(CGAL::rank_33<FT>(a2,b2,c2,a0,b0,c0,a1,b1,c1)==expected);
  assert(CGAL::rank_33<FT>(a2,b2,c2,a1,b1,c1,a0,b0,c0)==expected);
}

void test_rank_34(FT a0, FT b0, FT c0, FT d0,
                  FT a1, FT b1, FT c1, FT d1,
                  FT a2, FT b2, FT c2, FT d2,
                  int expected)
{
  std::cout << "testing:\n"
            << "\t " << a0 << "\t" << b0 << "\t" << c0 << "\t" << d0 << "\n"
            << "\t " << a1 << "\t" << b1 << "\t" << c1 << "\t" << d1 << "\n"
            << "\t " << a2 << "\t" << b2 << "\t" << c2 << "\t" << d2 << "\n";

  assert(CGAL::rank_34<FT>(a0,b0,c0,d0,a1,b1,c1,d1,a2,b2,c2,d2)==expected);
  assert(CGAL::rank_34<FT>(a0,b0,c0,d0,a2,b2,c2,d2,a1,b1,c1,d1)==expected);
  assert(CGAL::rank_34<FT>(a1,b1,c1,d1,a0,b0,c0,d0,a2,b2,c2,d2)==expected);
  assert(CGAL::rank_34<FT>(a1,b1,c1,d1,a2,b2,c2,d2,a0,b0,c0,d0)==expected);
  assert(CGAL::rank_34<FT>(a2,b2,c2,d2,a0,b0,c0,d0,a1,b1,c1,d1)==expected);
  assert(CGAL::rank_34<FT>(a2,b2,c2,d2,a1,b1,c1,d1,a0,b0,c0,d0)==expected);
}



int main()
{
  test_rank_33(1,0,0,
               0,1,0,
               0,0,1, 3);
  test_rank_33(1,0,0,
               0,1,0,
               1,0,1, 3);
  test_rank_33(1,0,0,
               1,1,1,
               1,0,1, 3);
  test_rank_33(1,0,0,
               1,1,1,
               1,0,1, 3);
  test_rank_33(1,0,0,
               1,1,0,
               1,5,0, 2);
  test_rank_33(0,0,0,
               0,0,0,
               0,0,0, 0);
  test_rank_33(1,2,3,
               1,2,3,
               1,2,3, 1);
  test_rank_33(1,2,3,
               2,4,6,
               4,8,12, 1);
  test_rank_33(1,2,3,
               2,4,6,
               4,8,11, 2);
  test_rank_33(1,0,1,
               1,0,1,
               1,0,1, 1);
  test_rank_33(1,1,1,
               1,0,1,
               1,0,1, 2);

  test_rank_33(0,1,1,
               0,0,1,
               0,0,1, 2);

  test_rank_34(1,0,0,0,
               0,1,0,0,
               0,0,1,0, 3);

  test_rank_34(1,0,0,0,
               0,1,0,0,
               1,1,0,0, 2);

  test_rank_34(1,0,0,0,
               0,1,0,0,
               1,1,0,1, 3);

  test_rank_34(0,1,0,0,
               0,0,1,0,
               0,1,1,0, 2);

  test_rank_34(0,1,0,0,
               0,0,1,0,
               0,1,1,1, 3);

  test_rank_34(1,0,0,0,
               0,0,1,0,
               1,0,1,0, 2);

  test_rank_34(1,0,0,0,
               0,0,1,0,
               1,0,1,1, 3);

  return 0;
}
