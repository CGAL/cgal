#include <CGAL/Cartesian.h>
#include <list>
#include <CGAL/Point_set_2.h>

typedef CGAL::Cartesian<double>     K;

typedef CGAL::Point_set_2<K>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2<K>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_2<K>                         Point;

CGAL::Point_set_2<K> PSet;
Point ar1[5];

int check_nn(const std::list<Vertex_handle>& L)
{
  std::cout << "check result!\n";
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
  std::list<Point> Lr;
  
  int w1,w2;
  Point p1(12,14);
  Point p2(-12,14);  
  Point p3(2,11);
  Point p4(5,6);
  Point p5(6.7,3.8);
  Point p6(11,20);
  Point p7(-5,6);  
  Point p8(12,0);
  Point p9(4,31);
  Point p10(-10,-10); 
  
  Lr.push_back(p1); Lr.push_back(p2); Lr.push_back(p3);
  Lr.push_back(p4); Lr.push_back(p5); Lr.push_back(p6);
  Lr.push_back(p7); Lr.push_back(p8); Lr.push_back(p9);
  Lr.push_back(p10); 
  
  PSet.insert(Lr.begin(),Lr.end()); 
  
   // init 
  ar1[0]=p4; ar1[1]=p5; ar1[2]=p3; ar1[3]=p7; ar1[4]=p8; 

  Point actual(30,45,10);

  // nearest neighbor ...  
  Vertex_handle v = PSet.nearest_neighbor(actual);
  std::cout << "Nearest neighbor:" << v->point() << "\n";
  
  if (v->point() == p4) w1=0; else w1=1;
  
  // k nearest neighbors ...
  std::list<Vertex_handle> L;
  std::list<Vertex_handle>::const_iterator it;

  PSet.nearest_neighbors(actual,5, std::back_inserter(L));
  std::cout << "actual point: " << actual << "\n";
    
  for (it=L.begin();it != L.end(); it++)
      std::cout << (*it)->point() << "\n";
     
   w2=check_nn(L);

  if (w1==0 && w2==0) return 0;
  else return 1;
}

