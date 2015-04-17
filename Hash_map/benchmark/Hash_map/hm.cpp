
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Timer.h>

#include <boost/unordered_map.hpp>


typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point_3;

typedef CGAL::Surface_mesh<Point_3>     Surface_mesh;
typedef CGAL::Polyhedron_3<Kernel>      Polyhedron_3;
typedef CGAL::Timer                     Timer;

template <typename G>
void
run(char* fname)
{
  typedef boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
  G g;
  std::ifstream input(fname);
  input >> g;

  typedef boost::property_map<G,CGAL::vertex_point_t>::type  VPM;
  VPM vpm = get(CGAL::vertex_point,g);


  std::vector<vertex_descriptor> V;
  std::vector<Point_3> P1, P2;
  BOOST_FOREACH(vertex_descriptor vd, vertices(g)){
    V.push_back(vd);
  }
  std::random_shuffle(V.begin(), V.end());

  std::map<vertex_descriptor,Point_3> vm;
  boost::unordered_map<vertex_descriptor,Point_3> vum;
  BOOST_FOREACH(vertex_descriptor vd, V){
    vm[vd] = get(vpm,vd);
    vum[vd] = get(vpm,vd);

  }

  Timer t;
  t.start();
  for(int i= 0; i < 1; i++){
  BOOST_FOREACH(vertex_descriptor vd, vertices(g)){
    P1.push_back(vm[vd]);
    /*
      if(vm[vd] != get(vpm,vd)){
      std::cerr << "error1" << std::endl;
    }
    */
  }
  }
  t.stop();
  std::cerr << t.time() << " sec."<< std::endl; 
  t.reset();
  t.start();
  for(int i= 0; i < 1; i++){
  BOOST_FOREACH(vertex_descriptor vd, vertices(g)){
    boost::unordered_map<vertex_descriptor,Point_3>::iterator it = vum.find(vd);
    P2.push_back((*it).second);
    /*
    if(vum[vd] != get(vpm,vd)){
      std::cerr << "error2" << std::endl;
    }
    */
  }
  }
  t.stop();
  std::cerr << t.time() << " sec."<< std::endl;
  if(P1.size() != P2.size()){
    std::cerr << "error3" << std::endl;
  } 
}

int main(int , char* argv[])
{
  run<Surface_mesh>(argv[1]);
  run<Polyhedron_3>(argv[1]);
  return 0;
}
