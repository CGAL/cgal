#include <list>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_2.h>

using namespace CGAL;
using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;

typedef CGAL::Point_set_2<K>::Edge           Edge;
typedef CGAL::Point_set_2<K>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2<K>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2<K>::Vertex         Vertex;

Point_set_2<K> PSet;

Point_2<K> ar1[5];

int check1(std::list<Vertex_handle> L)
{
  cout << "check 1!\n";
  if (L.size() != 5) return 1;
  std::list<Vertex_handle>::const_iterator it = L.begin();
  int i=0;
  int w=0;

  for(; it != L.end();it++){
    if (ar1[i] != (*it)->point()) w=1;
    i++;
  }
  return w;
}

int main()
{
  std::list<Point_2<K> > Lr;

  int w1,w2;

  Point_2<K> p1(12,14);
  Point_2<K> p2(-12,14);
  Point_2<K> p3(2,11);
  Point_2<K> p4(5,6);
  Point_2<K> p5(6.7,3.8);
  Point_2<K> p6(11,20);
  Point_2<K> p7(-5,6);
  Point_2<K> p8(12,0);
  Point_2<K> p9(4,31);
  Point_2<K> p10(-10,-10);

  Lr.push_back(p1); Lr.push_back(p2); Lr.push_back(p3);
  Lr.push_back(p4); Lr.push_back(p5); Lr.push_back(p6);
  Lr.push_back(p7); Lr.push_back(p8); Lr.push_back(p9);
  Lr.push_back(p10);

  PSet.insert(Lr.begin(),Lr.end());

   // init
  ar1[0]=p4; ar1[1]=p5; ar1[2]=p3; ar1[3]=p7; ar1[4]=p8;

  Point_2<K> actual(30,45,10);

  // nearest neighbor ...
  Vertex_handle v = PSet.nearest_neighbor(actual);
  std::cout << "Nearest neighbor:" << v->point() << "\n";

  if (v->point() == p4) w1=0; else w1=1;

  // k nearest neighbors ...
  std::list<Vertex_handle> L;
  std::list<Vertex_handle>::const_iterator it;

  PSet.nearest_neighbors(actual,5,std::back_inserter(L));
  std::cout << "actual point: " << actual << "\n";

  for (it=L.begin();it != L.end(); it++)
      std::cout << (*it)->point() << "\n";

   w2=check1(L);

   // construction ...
   Point_set_2<K> PSet2;
   Point_set_2<K> PSet3(Lr.begin(),Lr.end());
   Point_set_2<K> PSet4(Lr.begin(),Lr.end());

   // init ...
   PSet2.insert(Lr.begin(),Lr.end());

   PSet4.lookup(actual);

  if (w1==0 && w2==0) return 0;
  else return 1;
}
