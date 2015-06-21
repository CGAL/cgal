/*! \file test_connect_holes.cpp
 * Connecting a polygon with holes.
 */
 
#ifndef CGAL_BSO_RATIONAL_NT_H
#define CGAL_BSO_RATIONAL_NT_H

#include <CGAL/Exact_rational.h>
// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational         Number_type;
#endif

#include <CGAL/Cartesian.h>
#include <CGAL/connect_holes.h>
#include <list>
#include <iostream>
#include <string>
#include <sstream>

typedef CGAL::Cartesian<Number_type>               Kernel;
typedef Kernel::Point_2                            Point_2;
typedef CGAL::Polygon_2<Kernel>                    Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>         Polygon_with_holes_2;
typedef Polygon_2::FT										FT;
typedef Polygon_with_holes_2::Hole_const_iterator  Hole_const_iterator;
typedef Polygon_with_holes_2::Hole_iterator  		Hole_iterator;

//compute the area of a polygon with holes 
FT pwh_area(Polygon_with_holes_2 pwh)
{
  Polygon_2 outerP = pwh.outer_boundary();
  FT result = outerP.area();	
  if (! pwh.has_holes())
    return result;
  Hole_const_iterator hit = pwh.holes_begin();
  while (hit != pwh.holes_end()) {
    FT curHoleArea= (*hit).area();
    result = result + curHoleArea;
    ++hit;
  }
  //std::cout<< "The input pwh area is: " << result << std::endl;
  return result;
}

bool testExampleFile(const char* filename)
{  
  // Read a polygon with holes from a file.
  std::ifstream input_file (filename);
  if (! input_file.is_open())
  {
    std::cerr << "Failed to open the " << filename << std::endl;
    return false;
  }  
  Polygon_2               outerP;
  unsigned int            num_holes;

  input_file >> outerP;
  input_file >> num_holes;

  std::vector<Polygon_2>  holes (num_holes);
  unsigned int            k;

  for (k = 0; k < num_holes; k++)
    input_file >> holes[k];

  Polygon_with_holes_2    P (outerP, holes.begin(), holes.end());
  FT inputArea = pwh_area(P);
  
  // Connect the outer boundary of the polygon with its holes.
  std::list<Point_2>            pts;
  std::list<Point_2>::iterator  pit;

  connect_holes (P, std::back_inserter (pts));
  /*for (pit = pts.begin(); pit != pts.end(); ++pit)
    std::cout << '(' << *pit << ")  ";
  std::cout << std::endl;*/
  Polygon_2 res_polygon(pts.begin(),pts.end());
  FT outputArea = res_polygon.area();
  //std::cout<< "The output polygon area is: " << outputArea << std::endl;
  return (inputArea==outputArea);
}
/*
Test all 5 test files. The results are compared according to the area
signature. The input polygon with holes area is calculated, and then 
after the holes are connected, an output polygon is created, and its
area is calculated 
*/
int main()
{ 
  std::string testfilePrefix = "data/pgn_holes";
  std::string testfileSuffix = ".dat";   
  int result = 0;  
  for (int i = 1; i < 6; ++i) {
    std::stringstream strs;
    std::string si;
    strs << i;
    strs >> si;     
    std::string filename = testfilePrefix + si + testfileSuffix;
    const char * cfilename = filename.c_str();    
    bool res = testExampleFile(cfilename);
    if (!res) {
        std::cout << "test " << i << " was a bitter failure" << std::endl;
        result = 1;
    }  
    else {
      std::cout <<"test " << i << " was a great success" << std::endl;      
    }
  }
  if (result == 0)
    std::cout << "ALL TESTS SUCCEEDED!" << std::endl;  
  return result;
}

