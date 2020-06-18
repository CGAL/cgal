// test_output_modules.cpp

//----------------------------------------------------------
// Test the cgal environment for Optimal_transportation_reconstruction_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>

#include <iostream>
#include <cassert>
#include <vector>

#include "testing_tools.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_2                                          Point;
typedef K::Segment_2                                        Segment;

typedef CGAL::Optimal_transportation_reconstruction_2<K>    Otr_2;
typedef Otr_2::Vertex                                       Vertex;
typedef Otr_2::Rec_edge_2                                   R_edge_2;

void test_list_output(Otr_2& otr2);
void test_index_output(Otr_2& otr2);

int main ()
{
  std::vector<Point> points;
  //use the stair example for testing
  load_xy_file_points<Point>("data/stair-noise00.xy", points);

  Otr_2 otr2(points);
  otr2.run(100); //100 steps

  test_list_output(otr2);
  test_index_output(otr2);
}


void test_index_output(Otr_2& otr2)
{
  std::cout <<"(-------------OFF OUTPUT---------- )" << std::endl;

  std::vector<Point> points;
  std::vector<std::size_t> isolated_points;
  std::vector<std::pair<std::size_t,std::size_t> > edges;
  otr2.indexed_output(
      std::back_inserter(points),
      std::back_inserter(isolated_points),
      std::back_inserter(edges));

  std::stringstream sstr;

  sstr << "OFF " << points.size() <<
      " 0 " << edges.size()  << std::endl;

  for (std::vector<Point>::iterator it = points.begin();
       it != points.end(); it++)
  {
    sstr << *it << std::endl;
  }

  for (std::vector<std::pair<std::size_t,std::size_t> >::iterator
       it = edges.begin() ; it != edges.end() ; it++)
  {
    sstr << "2 " << it->first << " " << it->second << std::endl;
  }

  std::cout << sstr.str() << std::endl;

  //test cardinalities

  std::vector<std::string> res;
  for(;;)
  {
    std::string line;
    std::getline(sstr, line);
    if (!sstr.good())
      break;
    res.push_back(line);
  }

  assert(res.size() >= 60 && res.size() <= 105);
  assert(points.size() >= 35 && points.size() <= 70);
  assert(edges.size() >= 15 && edges.size() <= 40);

  for (std::size_t i = points.size() + 1 ; i < res.size() ; i++)
  {
    assert(res[i].substr(0,2) == "2 ");
  }
}

void test_list_output(Otr_2& otr2)
{
  std::cout <<"(-------------List OUTPUT---------- )" << std::endl;

  std::vector<Point> isolated_points;
  std::vector<Segment> edges;

  otr2.list_output(
    std::back_inserter(isolated_points), std::back_inserter(edges));

  int vertex_count = 0;
  for (std::vector<Point>::iterator it = isolated_points.begin();
       it != isolated_points.end(); it++)
  {
    vertex_count++;
    std::cout << *it << std::endl;
  }
  assert(vertex_count >= 8 && vertex_count <= 25);

  int edge_count = 0;
  for (std::vector<Segment>::iterator it = edges.begin();
       it != edges.end(); it++)
  {
    std::cout << *it << std::endl;
    edge_count++;
  }
  assert(edge_count >= 15 && edge_count <= 37);
}
