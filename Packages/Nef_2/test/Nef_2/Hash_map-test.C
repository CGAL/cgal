#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <list>
#include <CGAL/Hash_map.h>
#include <CGAL/test_macros.h>

using namespace CGAL;
using namespace std;
typedef list<int>::iterator list_handle;

int main() 
{
  CGAL_TEST_START;
 {
  list<int> L;
  L.push_back(1);
  list_handle it1 = L.begin();
  Hash_map<list_handle,int> H1,H2(-1);
  H1[it1] = 2; 
  CGAL_TEST(H1[it1]==2);
  CGAL_TEST(H2[it1]==-1);
  H1.clear(); H2.clear(-2);
  H2[it1] = 2; 
  CGAL_TEST(H1[it1]==0);
  CGAL_TEST(H2[it1]==2);
  list_handle it2 = L.end();
  const Hash_map<list_handle,int>* pH = &H2;
  CGAL_TEST((*pH)[it2]==-2);
 }
#if 0
 {
  point p1(1,2), p2(3,4);
  Hash_map<point,int,Kernel_index> H1,H2(-1);
  H1[p1] = 2; 
  CGAL_TEST(H1[p1]==2);
  CGAL_TEST(H2[p1]==-1);
  H1.clear(); H2.clear(-2);
  H2[p1] = 2; 
  CGAL_TEST(H1[p1]==0);
  CGAL_TEST(H2[p1]==2);
  const Hash_map<point,int,Kernel_index>* pH = &H2;
  CGAL_TEST((*pH)[p2]==-2);
  //H1.statistics();  H2.statistics();
 }
#endif
  CGAL_TEST_END;
  return 0;
}


