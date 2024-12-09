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
struct Point_3 : Kernel::Point_3 {
  using Kernel::Point_3::operator=;
  Point_3() : Kernel::Point_3(0,0,0) {}
  Point_3(int i) : Kernel::Point_3(i,0,0) {}
};
typedef Kernel::Vector_3                Vector_3;


typedef CGAL::Polyhedron_3<Kernel>      Mesh;
typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type  VPM;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef std::vector<vertex_descriptor> Vertex_list;

typedef CGAL::Timer                     Timer;

template <typename Map>
auto reserve(Map& sm, typename Map::size_type ii) -> decltype(sm.reserve(ii), void()){
    sm.reserve(ii);
}


template <typename Map>
Point_3 lookup(Map& sm,vertex_descriptor vh){
    auto search = sm.find(vh);
    if(search != sm.end())
        return search->second;
    return Point_3();
}

template <typename X, typename Y>
Point_3 lookup(CGAL::Unique_hash_map<X,Y>& sm,vertex_descriptor vh){
    const CGAL::Unique_hash_map<X,Y>& const_sm = sm;
    return const_sm[vh];
}

template <typename Map>
int fct(int ii, int jj, const Vertex_list& V, const Vertex_list& V2, const VPM& vpm, const std::string& s)
{
  int x = 0;
  Timer construct, query, lookups;
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

  int y=0;
  lookups.start();
  for(int j=0; j <jj; j++){
    for(vertex_descriptor vh : V2){
      y+= lookup(sm,vh).x();
    }
  }
  lookups.stop();
  if(y != 0) { std::cout << y << " != 0" << std::endl;}


  std::cerr << s << construct.time() << " sec.\t| " << query.time() << " sec.\t| " << lookups.time() << " sec." << std::endl;

  return x;
}

void random_mesh(int ii, int jj,Mesh& mesh,VPM& vpm,Vertex_list& V)
{

  for(int i =0; i < ii; i++){
    vertex_descriptor vd = add_vertex(mesh);
    put(vpm,vd, Point_3(i));
  }


  for(vertex_descriptor vd : vertices(mesh)){
    V.push_back(vd);
  }

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(V.begin(), V.end(), g);
}

void  fct(int ii, int jj)
{
  typedef std::map<vertex_descriptor,Point_3> SM;
  typedef std::unordered_map<vertex_descriptor,Point_3> SUM;
  typedef boost::unordered_map<vertex_descriptor,Point_3> BUM;
  typedef CGAL::Unique_hash_map<vertex_descriptor, Point_3> UHM;

  Mesh mesh1;
  VPM vpm1 = get(CGAL::vertex_point,mesh1);
  Vertex_list V1;
  random_mesh(ii,jj,mesh1,vpm1,V1);

  Mesh mesh2;
  VPM vpm2 = get(CGAL::vertex_point,mesh2);
  Vertex_list V2;
  random_mesh(ii,jj,mesh2,vpm2,V2);


  std::cerr << std::endl << ii << " items and queries (repeated " << jj << " times)" << std::endl;
  std::cerr << "Name\t\t\t| Version\t| Construction\t| Queries\t| Lookups" << std::endl;

  int temp;
  int res = fct<SM>(ii,jj, V1,V2, vpm1, "std::map\t\t|\t\t| " );
  temp = fct<SUM>(ii,jj,V1,V2, vpm1, "std::unordered_map\t|\t\t| " );
  if(temp != res){ std::cout << temp << " != " << res << std::endl;}
  temp = fct<BUM>(ii,jj, V1,V2, vpm1, "boost::unordered_map\t|\t\t| " );
  if(temp != res){ std::cout << temp << " != " << res << std::endl;}
  temp = fct<UHM>(ii,jj,V1,V2, vpm1, "CGAL::Unique_hash_map\t| master\t| " );
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
