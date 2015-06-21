#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <boost/foreach.hpp>


#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator   vertex_iterator;

typedef CGAL::Iterator_range<vertex_iterator> vertex_range;


vertex_range vertices_range(const Polyhedron& p)
{
  return vertex_range(vertices(p));
}

struct Fct
{
  void operator()(const vertex_descriptor& vd) const
  {
    std::cout << vd->point() << std::endl;
  }
};

void fct(const Polyhedron& p)
{
  vertex_range vr(vertices(p));
  
#ifndef CGAL_NO_CPP0X_RANGE_BASED_FOR
  std::cout << "new for loop" << std::endl;
  for(vertex_descriptor vd : vr){
    std::cout << vd->point() << std::endl;
  }
#endif
  
  std::cout << "BOOST_FOREACH" << std::endl;
  BOOST_FOREACH(vertex_descriptor vd, vr){
    std::cout << vd->point() << std::endl;
  }
  
  std::cout << "boost::tie + std::for_each" << std::endl;
  vertex_iterator vb, ve;
  
  boost::tie(vb,ve) = vertices_range(p);
  std::for_each(vb,ve, Fct());
}

int main(int, char** argv)
{
  Polyhedron P;  
  std::ifstream in(argv[1]);
  in >> P ;

  fct(P);
  return 0;
}
