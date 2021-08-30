#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Unique_hash_map.h>

#include <CGAL/Timer.h>
#include<boost/range/iterator_range.hpp>
#include <boost/unordered_map.hpp>
#include <unordered_map>

#include <boost/random/random_number_generator.hpp>
#include <boost/random/linear_congruential.hpp>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point_3;
typedef Kernel::Vector_3                Vector_3;


typedef CGAL::Polyhedron_3<Kernel>      Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

typedef CGAL::Timer                     Timer;



void  fct(int ii, int jj)
{
  typedef std::map<vertex_descriptor,Point_3> SM;
  typedef std::unordered_map<vertex_descriptor,Point_3> SUM;
  typedef boost::unordered_map<vertex_descriptor,Point_3> BUM;
  typedef CGAL::Unique_hash_map<vertex_descriptor, Point_3> UHM;


  Mesh mesh;
  typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type  VPM;
  VPM vpm = get(CGAL::vertex_point,mesh);


  Vector_3 v(0,0,0);

  for(int i =0; i < ii; i++){
    vertex_descriptor vd = add_vertex(mesh);
    put(vpm,vd, Point_3(i,0,0));
  }

  std::vector<vertex_descriptor> V;
  for(vertex_descriptor vd : vertices(mesh)){
    V.push_back(vd);
  }
  random_shuffle(V.begin(), V.end());


  Timer tsmc, tsumc, tbumc, tuhmc;
  Timer tsmq, tsumq, tbumq, tuhmq;

  for(int j=0; j <jj; j++){

    {
      tsmc.start();
      SM sm;
      for(vertex_descriptor vh : V){
        sm[vh] = get(vpm,vh);
      }
      tsmc.stop();


      tsmq.start();
      for(vertex_descriptor vh : V){
        v = v + (sm[vh] - CGAL::ORIGIN);
      }
      tsmq.stop();
    }
    {
      tsumc.start();
      SUM sm;
      for(vertex_descriptor vh : V){
        sm[vh] = get(vpm,vh);
      }
      tsumc.stop();


      tsumq.start();
      for(vertex_descriptor vh : V){
        v = v + (sm[vh] - CGAL::ORIGIN);
      }
      tsumq.stop();
    }
    {
      tbumc.start();
      BUM sm;
      for(vertex_descriptor vh : V){
        sm[vh] = get(vpm,vh);
      }
      tbumc.stop();


      tbumq.start();
      for(vertex_descriptor vh : V){
        v = v + (sm[vh] - CGAL::ORIGIN);
      }
      tbumq.stop();
    }

    {
      tuhmc.start();
      UHM sm;
      for(vertex_descriptor vh : V){
        sm[vh] = get(vpm,vh);
      }
      tuhmc.stop();


      tuhmq.start();
      for(vertex_descriptor vh : V){
        v = v + (sm[vh] - CGAL::ORIGIN);
      }
      tuhmq.stop();
    }
  }
  std::cerr << ii << " items and queries (repeated " << jj << " times)" << std::endl;
  std::cerr << "std::map             construction : "<< tsmc.time() << " sec." << std::endl;
  std::cerr << "std::map             queries      : "<< tsmq.time() << " sec." << std::endl;

  std::cerr << "std::unordered_map   construction : "<< tsumc.time() << " sec." << std::endl;
  std::cerr << "std::unordered_map   queries      : "<< tsumq.time() << " sec." << std::endl;

  std::cerr << "boost::unordered_map construction : "<< tbumc.time() << " sec." << std::endl;
  std::cerr << "boost::unordered_map queries      : "<< tbumq.time() << " sec.\n" << std::endl;

  std::cerr << "Unique_hash_map      construction : "<< tuhmc.time() << " sec." << std::endl;
  std::cerr << "Unique_hash_map      queries      : "<< tuhmq.time() << " sec.\n" << std::endl;

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
