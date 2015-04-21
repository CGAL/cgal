#define NOHASH 1

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
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Timer.h>
#include<boost/range/iterator_range.hpp> 
#include <boost/unordered_map.hpp>
#include <unordered_map>

#include <boost/random/random_number_generator.hpp>
#include <boost/random/linear_congruential.hpp>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point_3;
typedef Kernel::Vector_3                Vector_3;

typedef CGAL::Timer                     Timer;

template <typename G, typename Map>
void
run(const G& g)
{
  typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<G,CGAL::vertex_point_t>::const_type  VPM;
  VPM vpm = get(CGAL::vertex_point,g);

  std::vector<vertex_descriptor> V, V2;
  std::vector<Point_3> P1, P2;

  BOOST_FOREACH(vertex_descriptor vd, vertices(g)){
    V.push_back(vd);
    V2.push_back(vd);
  }
  
  boost::rand48 random;
  boost::random_number_generator<boost::rand48> rng(random);
  std::random_shuffle(V.begin(), V.end(), rng);
  
  Timer t;
#if 0
  t.start();
  Map vm;
  BOOST_FOREACH(vertex_descriptor vd, V){
    vm[vd] = get(vpm,vd);
  }
  t.stop();  std::cerr << "Insertion:  " << t.time() << " sec.     " << std::endl;
#endif

  Vector_3 v(0,0,0);
  std::size_t st=0;

#if 0
  std::cerr << "BOOST_FOREACH std::vector<vertex_descriptor)\n";
  t.reset(); t.start();
  for(int i=0; i<100; i++){
  BOOST_FOREACH(vertex_descriptor vd, V2){ 
#ifdef NOHASH    
    st += std::size_t(vd);
#else 
    typename Map::iterator it = vm.find(vd);
    v = v + ((*it).second - CGAL::ORIGIN);
#endif
  }
  }

  t.stop();  std::cerr << "  " <<t.time() << " sec.     " << std::endl;
#endif 

#if 1
  std::cerr << "BOOST_FOREACH boost::iterator_range r = vertices(g))\n";

  t.reset(); t.start();
  for(int i=0; i<100; i++){
  boost::iterator_range<typename boost::graph_traits<G>::vertex_iterator> r = vertices(g);

  BOOST_FOREACH(vertex_descriptor vd, r) {
#ifdef NOHASH    
    st += std::size_t(vd);
    // v = v + (get(vpm,vd) - CGAL::ORIGIN);
#else
    typename Map::iterator it = vm.find(vd);
    v = v + ((*it).second - CGAL::ORIGIN);
#endif
  }
  }
  t.stop();  std::cerr << "  " <<t.time() << " sec.     " << std::endl;
#endif


#if 1
   std::cerr << "BOOST_FOREACH CGAL::Iterator_range r = vertices(g))\n";
  t.reset(); t.start();
  for(int i=0; i<100; i++){
  CGAL::Iterator_range<typename boost::graph_traits<G>::vertex_iterator> ir = vertices(g);
  BOOST_FOREACH(vertex_descriptor vd, ir) {
#ifdef NOHASH
    st += std::size_t(vd);
    //v = v + (get(vpm,vd) - CGAL::ORIGIN);
#else
    typename Map::iterator it = vm.find(vd);
    v = v + ((*it).second - CGAL::ORIGIN);
#endif
  }
  }
  t.stop();  std::cerr << "  " <<t.time() << " sec.     " << std::endl;

#endif 


#if 1 
  std::cerr << "BOOST_FOREACH vertices(g))\n";
  t.reset(); t.start();
  for(int i=0; i<100; i++){
  BOOST_FOREACH(vertex_descriptor vd, vertices(g)) {
#ifdef NOHASH
    st += std::size_t(vd);
    //v = v + (get(vpm,vd) - CGAL::ORIGIN);
#else
    typename Map::iterator it = vm.find(vd);
    v = v + ((*it).second - CGAL::ORIGIN);
#endif
  }
  }
  t.stop();  std::cerr << "  " <<t.time() << " sec.     " << std::endl;
#endif 


#if 1
  std::cerr << "for vertices(g))\n";
  t.reset(); t.start();
  for(int i=0; i<100; i++){
    for(vertex_descriptor vd : vertices(g)) {
#ifdef NOHASH
    st += std::size_t(vd);
    //v = v + (get(vpm,vd) - CGAL::ORIGIN);
#else
    typename Map::iterator it = vm.find(vd);
    v = v + ((*it).second - CGAL::ORIGIN);
#endif
  }
  }
  t.stop();  std::cerr << "  " <<t.time() << " sec.     " << std::endl;
   
#endif

  
#if 0
  std::cerr << "boost::tie(vb,ve) = vertices(g);\n";
  t.reset(); t.start();
    for(int i=0; i<100; i++){
  typename boost::graph_traits<G>::vertex_iterator vb, ve;
  boost::tie(vb,ve) = vertices(g);
  for(; vb != ve; ++vb) {
    vertex_descriptor vd = *vb;
#ifdef NOHASH    
    st += std::size_t(vd);
    //v = v + (get(vpm,vd) - CGAL::ORIGIN);
#else
    typename Map::iterator it = vm.find(vd);
    v = v + ((*it).second - CGAL::ORIGIN);
#endif
  }
    }
  t.stop();  std::cerr << "  " <<t.time() << " sec.     " << std::endl;
#endif

  std::cerr << "v = " << v << "  " << st << std::endl;
  
}

struct blob {
  int a, b, c, d, e, f, g;
};


int main(int , char* argv[])
{

  {
    typedef CGAL::Surface_mesh<Point_3>     Mesh;
    typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
    typedef std::map<vertex_descriptor,Point_3> SM;
    typedef std::unordered_map<vertex_descriptor,Point_3> SUM;
    typedef boost::unordered_map<vertex_descriptor,Point_3> BUM;

    Mesh m;
    std::ifstream input(argv[1]);
    input >> m;

    std::cerr << num_vertices(m) << " items\n";
    std::cerr << "\nSurface_mesh  std::map"<< std::endl;
    run<Mesh,SM>(m);
#if 0
    std::cerr << "\nSurface_mesh  std::unordered_map"<< std::endl;
    run<Mesh,SUM>(m);
    std::cerr << "\nSurface_mesh  boost::unordered_map"<< std::endl;
    run<Mesh,BUM>(m);
#endif
  }

#if 0
  {
    typedef CGAL::Polyhedron_3<Kernel>      Mesh;
    typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
    typedef std::map<vertex_descriptor,Point_3> SM;
    typedef std::unordered_map<vertex_descriptor,Point_3> SUM;
    typedef boost::unordered_map<vertex_descriptor,Point_3> BUM;

    Mesh m;
    std::ifstream input(argv[1]);
    input >> m;

    std::cerr << "\nPolyhedron_3  std::map" << std::endl;
    run<Mesh,SM>(m);
    std::cerr << "\nPolyhedron_3 std::unordered_map"<< std::endl;
    run<Mesh,SUM>(m);
    std::cerr << "\nPolyhedron_3 boost::unordered_map"<< std::endl;
    run<Mesh,BUM>(m);    
  }
  
  {
    const int N = 3165798;
    std::cerr << "\nHashing "<< N << " pointers\n";
    std::vector<int> ints(N);

    std::vector<blob> data(N);
    for(int i =0; i <N ; i++){
      ints[i]=i;
    }
    std::random_shuffle(ints.begin(), ints.end());
    
    
    {
      boost::unordered_map<blob*,int> um;
      Timer t;
      t.start();
      for(int i= 0; i < N; i++){
        um[& (data[ints[i]])] = i;
      }
      t.stop();
      std::cerr << " boost::unordered_map: " <<  t.time() << " sec.\n";
    }
    
    {
      std::unordered_map<blob*,int> um;
      Timer t;
      t.start();
      for(int i= 0; i < N; i++){
        um[& (data[ints[i]])] = i;
      }
      t.stop();
      std::cerr  << " std::unordered_map: " << t.time() << " sec.\n";
    }
  }
#endif

  return 0;
}
