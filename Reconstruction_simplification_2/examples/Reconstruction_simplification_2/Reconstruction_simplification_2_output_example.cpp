#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>

#include<fstream>
#include<iostream>
#include <string>
#include <iterator>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::Segment_2                                        Segment;

typedef K::FT                                               FT;

typedef CGAL::Reconstruction_simplification_2<K> Rs_2;


void load_xy_file(const std::string& filename, std::list<Point>& points)
{
  std::ifstream ifs(filename);
  Point point;
  while (ifs >> point){
    points.push_back(point);
  }
  ifs.close();
}


void list_output(Rs_2& rs2)
{
  std::cout << "(-------------List output---------- )" << std::endl;
  
  std::vector<Point> isolated_vertices;
  std::vector<Segment> edges;
  
  rs2.extract_list_output(
    std::back_inserter(isolated_vertices), std::back_inserter(edges));
  
  for (std::vector<Point>::iterator vit = isolated_vertices.begin();
       vit != isolated_vertices.end(); vit++) {
    std::cout  <<  *vit << std::endl;
  }
  
  for (std::vector<Segment>::iterator eit = edges.begin();
       eit != edges.end(); eit++) {
    std::cout << *eit << std::endl;
  }
}


void index_output(Rs_2& rs2)
{  
  std::cout << "(-------------Off output---------- )" << std::endl;
  
  std::vector<Point> points;
  std::vector<std::size_t> isolated_points;
  std::vector<std::pair<std::size_t,std::size_t> > edges;

  rs2.indexed_output(
    std::back_inserter(points), 
    std::back_inserter(isolated_points),
    std::back_inserter(edges));
  
  std::cout << "OFF " << points.size() <<
		  " 0 " << edges.size()  << std::endl;

  for (std::vector<Point>::iterator it = points.begin(); 
       it != points.end(); it++) {
	  std::cout << *it << std::endl;
  }

  for (std::vector<std::pair<std::size_t,std::size_t> >::iterator it 
        = edges.begin();
		  it != edges.end(); it++) {
    std::cout << "2 "  << it->first << " " << it->second << std::endl;
  }
}

int main ()
{
  
  std::list<Point> points;
  
  load_xy_file("data/stair-noise00.xy", points);
  
  Rs_2 rs2(points);

  rs2.run(100); // 100 steps

  list_output(rs2);
  index_output(rs2);
  return 0;
}
