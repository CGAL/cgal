#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 400)
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
#include <CGAL/leda_rational.h>

using namespace CGAL;

typedef Cartesian<leda_rational> REP;
typedef point_set_traits_2<REP> TRAITS;
typedef Point_set_2<TRAITS>::Edge    Edge;
typedef Point_set_2<TRAITS>::Vertex  Vertex;

Point_set_2<TRAITS> PSet;

Point_2<REP> ar1[5];

int check1(std::list<Vertex> L)
{
  std::cout << "check 1!\n";
  if (L.size() != 5) return 1;
  std::list<Vertex>::const_iterator it = L.begin();
  int i=0;
  int w=0;
  
  for(; it != L.end();it++){
    if (ar1[i] != PSet.pos(*it)) w=1;
    i++;
  }
  return w;
}

int main()
{
  std::list<Point_2<REP> > Lr;
  
  int w1,w2;
  
  Point_2<REP> p1(12,14,1);
  Point_2<REP> p2(-12,14,1);  
  Point_2<REP> p3(2,11,1);
  Point_2<REP> p4(5,6,1);
  Point_2<REP> p5(67,38,10);
  Point_2<REP> p6(11,20,1);
  Point_2<REP> p7(-5,6,1);  
  Point_2<REP> p8(12,0,1);
  Point_2<REP> p9(4,31,1);
  Point_2<REP> p10(-10,-10,1); 
  
  Lr.push_back(p1); Lr.push_back(p2); Lr.push_back(p3);
  Lr.push_back(p4); Lr.push_back(p5); Lr.push_back(p6);
  Lr.push_back(p7); Lr.push_back(p8); Lr.push_back(p9);
  Lr.push_back(p10); 
  
  PSet.init(Lr.begin(),Lr.end()); 
  
   // init 
  ar1[0]=p4; ar1[1]=p5; ar1[2]=p3; ar1[3]=p7; ar1[4]=p8; 

  Point_2<REP> actual(30,45,10);

  // nearest neighbor ...  
  Vertex v = PSet.nearest_neighbor(actual);
  std::cout << "Nearest neighbor:" << PSet.pos(v) << "\n";
  
  if (PSet.pos(v) == p4) w1=0; else w1=1;
  
  // k nearest neighbors ...
  std::list<Vertex> L;
  std::list<Vertex>::const_iterator it;

  PSet.nearest_neighbors(actual,5,std::back_inserter(L));
  std::cout << "actual point: " << actual << "\n";
    
  for (it=L.begin();it != L.end(); it++)
      std::cout << PSet.pos(*it) << "\n";
     
   w2=check1(L);
   
   // construction ...
   Point_set_2<TRAITS> PSet2;
   Point_set_2<TRAITS> PSet3(Lr.begin(),Lr.end());
   Point_set_2<TRAITS> PSet4(Lr);
   
   // init ...
   PSet2.init(Lr);
   
   //get points/segments ...
   std::list<Point_2<REP> >    ptlist;
   std::list<Segment_2<REP> >  seglist;
   
   PSet3.points(std::back_inserter(ptlist));
   PSet3.segments(std::back_inserter(seglist));
   
   std::cout << PSet3.is_empty() << "\n";
   PSet3.clear();
   
   Edge eakt = PSet4.locate(actual);
   
   std::cout << PSet4.seg(eakt) << "\n";
   
   PSet4.lookup(actual);
   
   // insert/ delete
   Vertex vnew = PSet4.insert(actual);
   
   std::cout << PSet4.pos(vnew) << "\n";
   
   PSet4.del(actual);
 
   Edge he = PSet4.get_hull_edge();
   
   std::cout << PSet4.is_diagram_edge(he) << "\n";
   std::cout << PSet4.is_hull_edge(he) << "\n";
   
   std::cout << PSet4.dim() << "\n";
   
  if (w1==0 && w2==0) return 0;
  else return 1;
}

#endif
