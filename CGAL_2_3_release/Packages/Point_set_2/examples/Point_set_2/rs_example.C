#include <CGAL/Cartesian.h>
#include <list>
#include <CGAL/Point_set_2.h>

typedef CGAL::Cartesian<double>     K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Triangulation_face_base_2<K>   Fb;
typedef CGAL::Triangulation_default_data_structure_2<K,Vb,Fb> Tds;
typedef CGAL::Point_set_2<K,Tds>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_2<K>                         Point;

int main()
{
  CGAL::Point_set_2<K,Tds> PSet;
  std::list<Point> Lr;
  
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

  std::cout << "circular range search !\n";  
  CGAL::Circle_2<K> rc(p5,p6);

  std::list<Vertex_handle> LV;
  PSet.range_search(rc, std::back_inserter(LV));

  std::list<Vertex_handle>::const_iterator it;
  for (it=LV.begin();it != LV.end(); it++)
     std::cout << (*it)->point() << "\n";      
 
  std::cout << "triangular range search !\n";    
  
  LV.clear();
  PSet.range_search(p1,p2,p3, std::back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
     std::cout << (*it)->point() << "\n";    
  LV.clear();
 
  std::cout << "isorectangular range search !\n";
  Point pt1=p10; 
  Point pt3=p3; 
  Point pt2 = Point(pt3.x(),pt1.y());
  Point pt4 = Point(pt1.x(),pt3.y());
  
  PSet.range_search(pt1,pt2,pt3,pt4, std::back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
    std::cout << (*it)->point() << "\n"; 
  return 0;
}
