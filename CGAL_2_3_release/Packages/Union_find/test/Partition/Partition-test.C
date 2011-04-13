#include <cassert>
#include <CGAL/Partition.h>
#include <CGAL/test_macros.h>

typedef CGAL::Partition<int> Partition;
typedef Partition::item p_item;

int main() {
  CGAL_TEST_START;
  Partition P;
  p_item p1 = P.make_block(11);
  p_item p2 = P.make_block(-23);
  p_item p3 = P.make_block(17);
  P.union_blocks(p1,p2);
  CGAL_TEST(P.same_block(p1,p2));
  CGAL_TEST(!P.same_block(p1,p3));
  CGAL_TEST(P.size(p3)==1);
  CGAL_TEST(P.size(p1)==2);
  CGAL_TEST(P.inf(p1)==11);
  P.change_inf(p1,111);
  CGAL_TEST(P.inf(p1)==111);
  CGAL_TEST_END;
}

