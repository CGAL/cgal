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

  Hash_map<list_handle,int>::element_iterator it;
  for (it = H2.begin(); it != H2.end(); ++it)
  { *it; };
  CGAL_TEST_END;
}


