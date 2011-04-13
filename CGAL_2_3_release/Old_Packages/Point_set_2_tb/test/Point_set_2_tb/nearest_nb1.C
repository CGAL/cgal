#include <CGAL/basic.h>

#include <list>
#include <vector>
#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>

#include <CGAL/Point_set_2_tb.h>

using namespace CGAL;
using namespace std;

typedef Cartesian<double>     REP;
typedef CGAL::Triangulation_euclidean_traits_2<REP> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;


typedef point_set_traits_2<REP>                             TRAITS;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge           Edge;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex         Vertex;

Point_set_2_tb<TRAITS,Gt,Tds> PSet;

Point_2<REP> ar1[5];

int check1(std::list<Vertex_handle> L)
{
  cout << "check 1!\n";
  if (L.size() != 5) return 1;
  std::list<Vertex_handle>::const_iterator it = L.begin();
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
  
  Point_2<REP> p1(12,14);
  Point_2<REP> p2(-12,14);  
  Point_2<REP> p3(2,11);
  Point_2<REP> p4(5,6);
  Point_2<REP> p5(6.7,3.8);
  Point_2<REP> p6(11,20);
  Point_2<REP> p7(-5,6);  
  Point_2<REP> p8(12,0);
  Point_2<REP> p9(4,31);
  Point_2<REP> p10(-10,-10); 
  
  Lr.push_back(p1); Lr.push_back(p2); Lr.push_back(p3);
  Lr.push_back(p4); Lr.push_back(p5); Lr.push_back(p6);
  Lr.push_back(p7); Lr.push_back(p8); Lr.push_back(p9);
  Lr.push_back(p10); 
  
  PSet.init(Lr.begin(),Lr.end()); 
  
   // init 
  ar1[0]=p4; ar1[1]=p5; ar1[2]=p3; ar1[3]=p7; ar1[4]=p8; 

  Point_2<REP> actual(30,45,10);

  // nearest neighbor ...  
  Vertex_handle v = PSet.nearest_neighbor(actual);
  std::cout << "Nearest neighbor:" << PSet.pos(v) << "\n";
  
  if (PSet.pos(v) == p4) w1=0; else w1=1;
  
  // k nearest neighbors ...
  std::list<Vertex_handle> L;
  std::list<Vertex_handle>::const_iterator it;

  PSet.nearest_neighbors(actual,5,std::back_inserter(L));
  std::cout << "actual point: " << actual << "\n";
    
  for (it=L.begin();it != L.end(); it++)
      std::cout << PSet.pos(*it) << "\n";
     
   w2=check1(L);
   
   // construction ...
   Point_set_2_tb<TRAITS,Gt,Tds> PSet2;
   Point_set_2_tb<TRAITS,Gt,Tds> PSet3(Lr.begin(),Lr.end());
   Point_set_2_tb<TRAITS,Gt,Tds> PSet4(Lr);
   
   // init ...
   PSet2.init(Lr);
   
   //get points/segments ...
   std::list<Point_2<REP> >    ptlist;
   std::list<Segment_2<REP> >  seglist;
   
   PSet3.points(std::back_inserter(ptlist));
   PSet3.segments(std::back_inserter(seglist));
   
   std::cout << PSet3.is_empty() << "\n";
   PSet3.clear();
   
   PSet4.lookup(actual);
   
  if (w1==0 && w2==0) return 0;
  else return 1;
}
