
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>

#include <CGAL/Timer.h>
#include<boost/range/iterator_range.hpp>
#include <boost/unordered_map.hpp>
#include <unordered_map>

#include <boost/random/random_number_generator.hpp>
#include <boost/random/linear_congruential.hpp>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_2                 Point_2;
typedef Kernel::Vector_2                Vector_2;
typedef CGAL::Creator_uniform_2<double,Point_2>  Pt_creator;
typedef CGAL::Random_points_in_disc_2<Point_2,Pt_creator>  Random_points;;
typedef CGAL::Delaunay_triangulation_2<Kernel> Dt;
typedef Dt::Vertex_handle Vertex_handle;
typedef Dt::Finite_vertices_iterator Finite_vertices_iterator;

typedef CGAL::Timer                     Timer;





void  fct(int ii, int jj)
{
  typedef std::map<Vertex_handle,Point_2> SM;
  typedef std::unordered_map<Vertex_handle,Point_2> SUM;
  typedef boost::unordered_map<Vertex_handle,Point_2> BUM;

  Dt dt;
  Vector_2 v(0,0);
  Random_points rp( 250);
  std::vector<Point_2> points;
  for(int i =0; i < ii; i++){
    Point_2 p = *rp++;
    points.push_back(p);
  }

  dt.insert(points.begin(), points.end());

  std::vector<Vertex_handle> vertices;
  Finite_vertices_iterator b = dt.finite_vertices_begin(), e = dt.finite_vertices_end();
  for(; b!=e; ++b){
    vertices.push_back(b);
  }
  random_shuffle(vertices.begin(), vertices.end());

  Timer tsmc, tsumc, tbumc;
  Timer tsmq, tsumq, tbumq;
  for(int j=0; j <jj; j++){

    {
      tsmc.start();
      SM sm;
      for(Vertex_handle vh : vertices){
        sm[vh] = vh->point();
      }
      tsmc.stop();


      tsmq.start();
      for(Vertex_handle vh : vertices){
        v = v + (sm[vh] - CGAL::ORIGIN);
      }
      tsmq.stop();
    }
    {
      tsumc.start();
      SUM sm;
      for(Vertex_handle vh : vertices){
        sm[vh] = vh->point();
      }
      tsumc.stop();


      tsumq.start();
      for(Vertex_handle vh : vertices){
        v = v + (sm[vh] - CGAL::ORIGIN);
      }
      tsumq.stop();
    }
    {
      tbumc.start();
      BUM sm;
      for(Vertex_handle vh : vertices){
        sm[vh] = vh->point();
      }
      tbumc.stop();


      tbumq.start();
      for(Vertex_handle vh : vertices){
        v = v + (sm[vh] - CGAL::ORIGIN);
      }
      tbumq.stop();
    }
  }
  std::cerr << ii << " items and queries (repeated " << jj << " times)" << std::endl;
  std::cerr << "std::map             construction : "<< tsmc.time() << " sec." << std::endl;
  std::cerr << "std::map             queries      : "<< tsmq.time() << " sec." << std::endl;

  std::cerr << "std::unordered_map   construction : "<< tsumc.time() << " sec." << std::endl;
  std::cerr << "std::unordered_map   queries      : "<< tsumq.time() << " sec." << std::endl;

  std::cerr << "boost::unordered_map construction : "<< tbumc.time() << " sec." << std::endl;
  std::cerr << "boost::unordered_map queries      : "<< tbumq.time() << " sec.\n" << std::endl;
}

int main(int , char* argv[])
{
  fct(1000000, 10);
  fct(100000, 100);
  fct(10000, 1000);
  fct(1000, 10000);
  fct(500,  20000);
  fct(250,  40000);
  fct(100, 100000);
  return 0;
}
