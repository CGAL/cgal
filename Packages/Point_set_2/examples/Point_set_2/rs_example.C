#include <CGAL/basic.h>
#include <list>
#include <vector>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_set_2.h>

using namespace CGAL;
using namespace std;

typedef Cartesian<double>     K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Triangulation_face_base_2<K>   Fb;
typedef CGAL::Triangulation_default_data_structure_2<K,Vb,Fb> Tds;
typedef CGAL::Point_set_2<K,Tds>::Vertex_handle  Vertex_handle;

int main()
{
  Point_set_2<K,Tds> PSet;
  std::list<Point_2<K> > Lr;
  
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

  PSet.init(Lr.begin(),Lr.end()); 

  cout << "circular range search !\n";  
  Circle_2<K> rc(p5,p6);

  std::list<Vertex_handle> LV;
  PSet.range_search(rc,back_inserter(LV));

  std::list<Vertex_handle>::const_iterator it;
  for (it=LV.begin();it != LV.end(); it++)
     cout << PSet.pos(*it) << "\n";      
 
  cout << "triangular range search !\n";    
  
  LV.clear();
  PSet.range_search(p1,p2,p3,back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
     cout << PSet.pos(*it) << "\n";    
  LV.clear();
 
  cout << "isorectangular range search !\n";
  Point_2<K> pt1=p10; 
  Point_2<K> pt3=p3; 
  Point_2<K> pt2 = Point_2<K>(pt3.x(),pt1.y());
  Point_2<K> pt4 = Point_2<K>(pt1.x(),pt3.y());
  
  PSet.range_search(pt1,pt2,pt3,pt4,back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
    cout << PSet.pos(*it) << "\n"; 
  return 0;
}
