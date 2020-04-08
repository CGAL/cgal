#include <CGAL/Simple_cartesian.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>


#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>


typedef CGAL::Simple_cartesian<double>              Kernel;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> LCC_traits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
         <2, 3, LCC_traits>::type LCC;

typedef boost::graph_traits<LCC>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<LCC>::vertex_iterator   vertex_iterator;

typedef CGAL::Iterator_range<vertex_iterator> vertex_range;

vertex_range vertices_range(const LCC& lcc)
{
  return vertex_range(vertices(lcc));
}

struct Fct
{
  void operator()(const vertex_descriptor& vd) const
  {
    std::cout << vd->point() << std::endl;
  }
};

void fct(const LCC& lcc)
{
  vertex_range vr(vertices(lcc));

  std::cout << "new for loop" << std::endl;
  for(vertex_descriptor vd : vr){
    std::cout << vd->point() << std::endl;
  }

  std::cout << "boost::tie + std::for_each" << std::endl;
  vertex_iterator vb, ve;

  boost::tie(vb,ve) = vertices_range(lcc);
  std::for_each(vb,ve, Fct());
}

int main(int argc, char** argv)
{
  LCC lcc;
  CGAL::read_off((argc>1)?argv[1]:"cube.off", lcc);

  fct(lcc);
  return 0;
}
