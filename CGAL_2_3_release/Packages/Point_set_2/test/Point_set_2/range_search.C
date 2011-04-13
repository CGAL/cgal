#include <CGAL/basic.h>
#include <list>
#include <vector>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_set_2.h>

using namespace CGAL;
using namespace std;

typedef Cartesian<double>     K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Triangulation_face_base_2<K>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<K,Vb,Fb> Tds;

typedef CGAL::Point_set_2<K,Tds>::Edge    Edge;
typedef CGAL::Point_set_2<K,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2<K,Tds>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2<K,Tds>::Vertex  Vertex;

Point_set_2<K,Tds> PSet;

Point_2<K> ar1[6];
Point_2<K> ar2[3];
Point_2<K> ar3[3];

int check1(std::list<Vertex_handle> L)
{
  cout << "check 1!\n";
  if (L.size() != 6) return 1;
  std::list<Vertex_handle>::const_iterator it = L.begin();
  int i=0;
  int w=0;
  
  for(; it != L.end();it++){
    if (ar1[i] != (*it)->point()) w=1;
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
    if (ar2[i] != (*it)->point()) w=1;
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
    if (ar3[i] != (*it)->point()) w=1;
    i++;
 }
 
  if (w==1) cout << "check3 failed!\n";
 return w;
}

int main()
{
  Point_2<K> pnew(12,6.2);
  
  int w1,w2,w3;

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
  
  // init 
  ar1[0]=p1; ar1[1]=p6; ar1[2]=p3; ar1[3]=p4; ar1[4]=p5; ar1[5]=pnew; 
  ar2[0]=p1; ar2[1]=p3; ar2[2]=p2;
  ar3[0]=p7; ar3[1]=p10; ar3[2]=p3;
  
  Lr.push_back(p1); Lr.push_back(p2); Lr.push_back(p3);
  Lr.push_back(p4); Lr.push_back(p5); Lr.push_back(p6);
  Lr.push_back(p7); Lr.push_back(p8); Lr.push_back(p9);
  Lr.push_back(p10); 

  PSet.insert(Lr.begin(),Lr.end()); 

  cout << "insert!\n"; cout.flush();
  PSet.insert(pnew);


  cout << "range search for circle !\n";  
  Circle_2<K> rc(p5,p6);

  std::list<Vertex_handle> LV;
  PSet.range_search(rc,back_inserter(LV));

  std::list<Vertex_handle>::const_iterator it;
  for (it=LV.begin();it != LV.end(); it++)
     cout << (*it)->point() << "\n";      
     
  w1 = check1(LV);
 
  cout << "range search for triangle !\n";    
  
  LV.clear();
  PSet.range_search(p1,p2,p3,back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
     cout << (*it)->point() << "\n";
    
  w2 = check2(LV);    
  LV.clear();
 
  cout << "range search for iso rectangle !\n";
  Point_2<K> pt1=p10; // lower left
  Point_2<K> pt3=p3; // upper right 
  
  Point_2<K> pt2 = Point_2<K>(pt3.x(),pt1.y());
  Point_2<K> pt4 = Point_2<K>(pt1.x(),pt3.y());
  
  PSet.range_search(pt1,pt2,pt3,pt4,back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
    cout << (*it)->point() << "\n"; 

  w3 = check3(LV);
  
  if (w1==0 && w2==0 && w3==0) return 0;
  else return 1;
}
