#include <cassert>
#include <CGAL/Union_find.h>
#include <CGAL/test_macros.h>

typedef CGAL::Union_find<int> Union_find;
typedef Union_find::handle handle;
typedef Union_find::iterator Iterator;
typedef Union_find::const_handle const_handle;
typedef Union_find::const_iterator const_iterator;

int main() {
  CGAL_TEST_START;
  Union_find P;
  handle p1 = P.make_set(11);
  handle p2 = P.make_set(-23);
  handle p3 = P.push_back(17);
  CGAL_TEST(P.size()==3);
  CGAL_TEST(P.number_of_sets()==3);
  P.unify_sets(p1,p2);
  CGAL_TEST(P.size()==3);
  CGAL_TEST(P.number_of_sets()==2);
  CGAL_TEST(P.same_set(p1,p2));
  CGAL_TEST(!P.same_set(p1,p3));
  CGAL_TEST(P.size(p3)==1);
  CGAL_TEST(P.size(p1)==2);
  CGAL_TEST(*p1 == 11);
  *p1 = 111;
  CGAL_TEST(*p1 == 111);
  CGAL_TEST(P.size()==3);
  CGAL_TEST(P.bytes()>0);
  int i = 0;
  for(Iterator it = P.begin(); it != P.end(); ++it) *it = i++;
  int A[] = {1,2,3};
  P.insert(A,A+3);
  CGAL_TEST_END;
}

