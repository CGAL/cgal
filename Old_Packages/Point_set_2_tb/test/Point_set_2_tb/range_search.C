
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


typedef point_set_traits_2<REP>      TRAITS;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge    Edge;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex  Vertex;


Point_set_2_tb<TRAITS,Gt,Tds> PSet;

Point_2<REP> ar1[6];
Point_2<REP> ar2[3];
Point_2<REP> ar3[3];

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
  
  if (w==1) cout << "check1 failed!\n";
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
  if (w==1) cout << "check2 failed!\n";
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
 
  if (w==1) cout << "check3 failed!\n";
 return w;
}

int main()
{
  Point_2<REP> pnew(12,6.2);
  
  int w1,w2,w3;

  std::list<Point_2<REP> > Lr;
  
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
  
  // init 
  ar1[0]=p1; ar1[1]=p6; ar1[2]=p3; ar1[3]=p4; ar1[4]=p5; ar1[5]=pnew; 
  ar2[0]=p1; ar2[1]=p3; ar2[2]=p2;
  ar3[0]=p7; ar3[1]=p10; ar3[2]=p3;
  
  Lr.push_back(p1); Lr.push_back(p2); Lr.push_back(p3);
  Lr.push_back(p4); Lr.push_back(p5); Lr.push_back(p6);
  Lr.push_back(p7); Lr.push_back(p8); Lr.push_back(p9);
  Lr.push_back(p10); 

  PSet.init(Lr.begin(),Lr.end()); 

  cout << "insert!\n"; cout.flush();
  PSet.insert(pnew);


  cout << "range search for circle !\n";  
  Circle_2<REP> rc(p5,p6);

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
  Point_2<REP> pt1=p10; // lower left
  Point_2<REP> pt3=p3; // upper right 
  
  Point_2<REP> pt2 = Point_2<REP>(pt3.x(),pt1.y());
  Point_2<REP> pt4 = Point_2<REP>(pt1.x(),pt3.y());
  
  PSet.range_search(pt1,pt2,pt3,pt4,back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
    cout << PSet.pos(*it) << "\n"; 

  w3 = check3(LV);
  
  if (w1==0 && w2==0 && w3==0) return 0;
  else return 1;
}
