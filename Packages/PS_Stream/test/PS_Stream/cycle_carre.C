#include <CGAL/Cartesian.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <CGAL/IO/Color.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/IO/PS_Stream_3.h>

typedef CGAL::Cartesian<double> D;
typedef CGAL::Bbox_3            PS_BBox3;
typedef D::Direction_3          Direction;
typedef D::Point_3              Point3;


int main()
{
  double x,y,z,lx,ly,lz;
  std::string filename;
  std::cout << "Enter file name: ";std::cin >> filename;
  std::cerr << "Enter view x: ";std::cin >> x;
  std::cerr << "Enter view y: ";std::cin >> y;
  std::cerr << "Enter view z: ";std::cin >> z;
  std::cerr << "Enter light x: ";std::cin >> lx;
  std::cerr << "Enter light y: ";std::cin >> ly;
  std::cerr << "Enter light z: ";std::cin >> lz;

  Direction dir(x,y,z);
  Direction light(lx,ly,lz);
  PS_BBox3 bb3(-6,-6,-6,6,6,6);
  CGAL::PS_Stream_3 ps(bb3,dir,light,300,filename.c_str(), CGAL::PS_Stream::READABLE_EPS);
  
  //Cycle
  Point3 a1(-3,-5,-1);Point3 b1(-3,5,1);Point3 c1(-4,5,1);Point3 d1(-4,-5,-1);
  
  vector<Point3> v1;  
  v1.push_back(a1);v1.push_back(b1);v1.push_back(c1);v1.push_back(d1);
  vector<Point3> v1bis;  
  v1bis.push_back(d1);v1bis.push_back(c1);v1bis.push_back(b1);v1bis.push_back(a1);
  
  Point3 a2(-5,3,0);Point3 b2(5,3,0);Point3 c2(5,4,0);Point3 d2(-5,4,0);
  vector<Point3> v2;
  v2.push_back(a2);v2.push_back(b2);v2.push_back(c2);v2.push_back(d2);
  vector<Point3> v2bis;  
  v2bis.push_back(d2);v2bis.push_back(c2);v2bis.push_back(b2);v2bis.push_back(a2);

  Point3 a3(4,-5,1);Point3 b3(4,5,-1);Point3 c3(3,5,-1);Point3 d3(3,-5,1);
  vector<Point3> v3;
  v3.push_back(a3);v3.push_back(b3);v3.push_back(c3);v3.push_back(d3);
  vector<Point3> v3bis;  
  v3bis.push_back(d3);v3bis.push_back(c3);v3bis.push_back(b3);v3bis.push_back(a3);

  Point3 a4(-5,-4,0);Point3 b4(5,-4,0);Point3 c4(5,-3,0);Point3 d4(-5,-3,0);
  vector<Point3> v4;
  v4.push_back(a4);v4.push_back(b4);v4.push_back(c4);v4.push_back(d4);
  vector<Point3> v4bis;  
  v4bis.push_back(d4);v4bis.push_back(c4);v4bis.push_back(b4);v4bis.push_back(a4);

  CGAL::PS_facet_3 face1(v1,CGAL::BLACK,CGAL::RED,CGAL::UNIFORM_FILL,1);
  CGAL::PS_facet_3 face2(v2,CGAL::BLACK,CGAL::GREEN,CGAL::UNIFORM_FILL,2);
  CGAL::PS_facet_3 face3(v3,CGAL::BLACK,CGAL::BLUE,CGAL::UNIFORM_FILL,3); 
  CGAL::PS_facet_3 face4(v4,CGAL::BLACK,CGAL::ORANGE,CGAL::UNIFORM_FILL,4); 
 
  CGAL::PS_facet_3 face1bis(v1bis,CGAL::BLACK,CGAL::BLACK,CGAL::UNIFORM_FILL,5);
  CGAL::PS_facet_3 face2bis(v2bis,CGAL::BLACK,CGAL::BLUE,CGAL::UNIFORM_FILL,6);
  CGAL::PS_facet_3 face3bis(v3bis,CGAL::BLACK,CGAL::VIOLET,CGAL::UNIFORM_FILL,7); 
  CGAL::PS_facet_3 face4bis(v4bis,CGAL::BLACK,CGAL::RED,CGAL::UNIFORM_FILL,8); 

  vector<CGAL::PS_facet_3> vfacet;
  vfacet.push_back(face1);vfacet.push_back(face2);vfacet.push_back(face3);vfacet.push_back(face4); 
  vfacet.push_back(face1bis);vfacet.push_back(face2bis);vfacet.push_back(face3bis);vfacet.push_back(face4bis); 
  
    for(unsigned int i=0;i<vfacet.size();i++) {
    ps.add_facet(vfacet[i]);
  }

  ps.display();
  return 0;
}
