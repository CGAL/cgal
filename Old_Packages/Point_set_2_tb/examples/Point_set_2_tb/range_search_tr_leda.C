
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
#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_leda_traits_2.h>

#include <CGAL/point_set_leda_traits_2.h>
#include <CGAL/Point_set_2_tb.h>
#include <CGAL/leda_rational.h>

using namespace CGAL;
using namespace std;

//typedef CGAL::Triangulation_euclidean_leda_rat_traits_2 Gt;
typedef CGAL::Triangulation_euclidean_leda_float_traits_2 Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;

//typedef CGAL::point_set_leda_rat_traits_2      TRAITS;
typedef CGAL::point_set_leda_float_traits_2      TRAITS;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge    Edge;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex  Vertex;
typedef Gt::Segment_2 Segment;
typedef Gt::Point_2 Point;
//typedef leda_rat_circle Circle;
typedef leda_circle Circle;


Point_set_2_tb<TRAITS,Gt,Tds> PSet;

Point ar1[6];
Point ar2[3];
Point ar3[3];

int check1(std::list<Vertex_handle> L)
{
  cout << "check 1!\n";
  if (L.size() != 6) return 1;
  std::list<Vertex_handle>::const_iterator it = L.begin();
  int i=0;
  int w=0;
  
  for(; it != L.end();it++){
    if (ar1[i] != PSet.pos(*it)) w=1;
    i++;
  }
  return w;
}

int check2(std::list<Vertex_handle> L)
{
  cout << "check 2!\n";
  if (L.size() != 3) return 1; 
  std::list<Vertex_handle>::const_iterator it = L.begin();
  int i=0;
  int w=0;
    
  for(; it != L.end();it++){
    if (ar2[i] != PSet.pos(*it)) w=1;
    i++;
  }
  return w;
}

int check3(std::list<Vertex_handle> L)
{
 cout << "check 3!\n";
 if (L.size() != 3) return 1; 
 std::list<Vertex_handle>::const_iterator it = L.begin();
 int i=0;
 int w=0;
    
 for(; it != L.end();it++){
    if (ar3[i] != PSet.pos(*it)) w=1;
    i++;
 }
 return w;
}

int main()
{
  Point pnew(120,62,10);
  
  int w1,w2,w3;

  std::list<Point>  Lr;
  
  Point p1(12,14,1);
  Point p2(-12,14,1);  
  Point p3(2,11,1);
  Point p4(5,6,1);
  Point p5(67,38,10);
  Point p6(11,20,1);
  Point p7(-5,6,1);  
  Point p8(12,0,1);
  Point p9(4,31,1);
  Point p10(-10,-10,1); 
  
  // init 
  ar1[0]=p1; ar1[1]=p6; ar1[2]=p3; ar1[3]=p4; ar1[4]=p5; ar1[5]=pnew; 
  ar2[0]=p2; ar2[1]=p3; ar2[2]=p1;
  ar3[0]=p7; ar3[1]=p10; ar3[2]=p3;
  
  Lr.push_back(p1); Lr.push_back(p2); Lr.push_back(p3);
  Lr.push_back(p4); Lr.push_back(p5); Lr.push_back(p6);
  Lr.push_back(p7); Lr.push_back(p8); Lr.push_back(p9);
  Lr.push_back(p10); 

  PSet.init(Lr.begin(),Lr.end()); 

  cout << "insert!\n"; cout.flush();
  PSet.insert(pnew);


  cout << "range search for circle !\n";  
  Circle rc(p5,p6);

  std::list<Vertex_handle> LV;
  PSet.range_search(rc,back_inserter(LV));

  std::list<Vertex_handle>::const_iterator it;
  for (it=LV.begin();it != LV.end(); it++)
     cout << PSet.pos(*it) << "\n";      
     
  w1 = check1(LV);
 
  cout << "range search for triangle !\n";    
  
  LV.clear();
  PSet.range_search(p1,p2,p3,back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
     cout << PSet.pos(*it) << "\n";
    
  w2 = check2(LV);    
  LV.clear();
 
  cout << "range search for iso rectangle !\n";
  Point pt1=p10; // lower left
  Point pt3=p3; // upper right 
  
  Point pt2 = Point(pt3.xcoord(),pt1.ycoord());
  Point pt4 = Point(pt1.xcoord(),pt3.ycoord());
  
  PSet.range_search(pt1,pt2,pt3,pt4,back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
    cout << PSet.pos(*it) << "\n"; 

  w3 = check3(LV);
  
  if (w1==0 && w2==0 && w3==0) return 0;
  else return 1;
}

#endif
