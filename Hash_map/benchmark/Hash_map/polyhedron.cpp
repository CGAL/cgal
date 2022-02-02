#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <random>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Unique_hash_map.h>

#include <CGAL/Timer.h>
#include<boost/range/iterator_range.hpp>
#include <boost/unordered_map.hpp>
#include <unordered_map>

#include <boost/random/random_number_generator.hpp>
#include <boost/random/linear_congruential.hpp>

typedef CGAL::Simple_cartesian<int>  Kernel;
typedef Kernel::Point_3                 Point_3;
typedef Kernel::Vector_3                Vector_3;


typedef CGAL::Polyhedron_3<Kernel>      Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

typedef CGAL::Timer                     Timer;

template <typename Map>
auto reserve(Map& sm, typename Map::size_type ii) -> decltype(sm.reserve(ii), void()){
    sm.reserve(ii);
}

template <typename Map, typename VPM>
double fct(int ii, int jj, const std::vector<vertex_descriptor>& V, const VPM& vpm, const std::string& s)
{
  int x = 0;
  Timer construct, query;
  construct.start();
  for(int j=0; j <jj; j++){
    Map sm;
    reserve(sm,ii);
    for(vertex_descriptor vh : V){
      sm[vh] = get(vpm,vh);
    }
  }
  construct.stop();

  Map sm;
  for(vertex_descriptor vh : V){
    sm[vh] = get(vpm,vh);
  }

  query.start();
  for(int j=0; j <jj; j++){
    for(vertex_descriptor vh : V){
      x+= sm[vh].x();
    }
  }
  query.stop();


  std::cerr << s << " construction : "<< construct.time() << " sec.    ";
  std::cerr      << " queries      : "<< query.time()     << " sec." << std::endl;

  return x;
}



void  fct(int ii, int jj)
{
  typedef std::map<vertex_descriptor,Point_3> SM;
  typedef std::unordered_map<vertex_descriptor,Point_3> SUM;
  typedef boost::unordered_map<vertex_descriptor,Point_3> BUM;
  typedef CGAL::Unique_hash_map<vertex_descriptor, Point_3> UHM;


  Mesh mesh;
  typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type  VPM;
  VPM vpm = get(CGAL::vertex_point,mesh);

  for(int i =0; i < ii; i++){
    vertex_descriptor vd = add_vertex(mesh);
    put(vpm,vd, Point_3(i,0,0));
  }

  std::vector<vertex_descriptor> V;
  for(vertex_descriptor vd : vertices(mesh)){
    V.push_back(vd);
  }

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(V.begin(), V.end(), g);


    std::cerr << std::endl << ii << " items and queries (repeated " << jj << " times)" << std::endl;

  int temp;
  int res = fct<SM>(ii,jj, V, vpm, "std::map             " );
  temp = fct<SUM>(ii,jj,V, vpm, "std::unordered_map   " );
  if(temp != res){ std::cout << temp << " != " << res << std::endl;}
  temp = fct<BUM>(ii,jj, V, vpm, "boost::unordered_map " );
  if(temp != res){ std::cout << temp << " != " << res << std::endl;}
  temp = fct<UHM>(ii,jj, V, vpm, "CGAL::Unique_hash_map" );
  if(temp != res){ std::cout << temp << " != " << res << std::endl;}
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
  fct(10, 1000000);
  return 0;
}
