#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA installed!\n";
 return 0;
}
#else 

#include <CGAL/config.h>
#include <list>
#include <vector>
#include <CGAL/Point_set_2.h>

typedef CGAL::Cartesian<double>          REP;
typedef CGAL::point_set_traits_2<REP>      TRAITS;
typedef CGAL::Point_set_2<TRAITS>::Edge    Edge;
typedef CGAL::Point_set_2<TRAITS>::Vertex  Vertex;

int main()
{
  CGAL::Point_set_2<TRAITS> PSet;
  CGAL::Point_2<REP> pnew;

  std::list<CGAL::Point_2<REP> > Lr;
  
  CGAL::Point_2<REP> p1(12,14);
  CGAL::Point_2<REP> p2(-12,14);  
  CGAL::Point_2<REP> p3(2,11);
  CGAL::Point_2<REP> p4(5,6);
  CGAL::Point_2<REP> p5(6.7,3.8);
  CGAL::Point_2<REP> p6(11,20);
  CGAL::Point_2<REP> p7(-5,6);  
  CGAL::Point_2<REP> p8(12,0);
  CGAL::Point_2<REP> p9(4,31);
  CGAL::Point_2<REP> p10(-10,-10); 
  
  Lr.push_back(p1); Lr.push_back(p2); Lr.push_back(p3);
  Lr.push_back(p4); Lr.push_back(p5); Lr.push_back(p6);
  Lr.push_back(p7); Lr.push_back(p8); Lr.push_back(p9);
  Lr.push_back(p10); 

  PSet.init(Lr.begin(),Lr.end()); 

  pnew= CGAL::Point_2<REP>(12,6.2);
  std::cout << "insert!\n"; 
  PSet.insert(pnew);

  std::cout << "range search for circle !\n";  
  CGAL::Circle_2<REP> rc(p5,p6);

  std::list<Vertex> LV;
  PSet.range_search(rc,std::back_inserter(LV));

  std::list<Vertex>::const_iterator it;
  for (it=LV.begin();it != LV.end(); it++)
     std::cout << PSet.pos(*it) << "\n";      
 
  std::cout << "range search for triangle !\n";    
  
  LV.clear();
  PSet.range_search(p1,p2,p3,std::back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
     std::cout << PSet.pos(*it) << "\n";
        
  LV.clear();
 
  std::cout << "range search for iso rectangle !\n";
  CGAL::Point_2<REP> pt1=p10; // lower left
  CGAL::Point_2<REP> pt3=p3; // upper right 
  
  CGAL::Point_2<REP> pt2 = CGAL::Point_2<REP>(pt3.x(),pt1.y());
  CGAL::Point_2<REP> pt4 = CGAL::Point_2<REP>(pt1.x(),pt3.y());
  
  PSet.range_search(pt1,pt2,pt3,pt4,std::back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
    std::cout << PSet.pos(*it) << "\n"; 

  return 0;
}
#endif
