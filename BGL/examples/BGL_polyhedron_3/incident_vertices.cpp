#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel>     Polyhedron;

typedef boost::graph_traits<Polyhedron> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;
typedef CGAL::Halfedge_around_target_iterator<Polyhedron> halfedge_around_target_iterator;

template <typename OutputIterator>
OutputIterator
adjacent_vertices_V1(const Polyhedron& g,
                     vertex_descriptor vd,
                     OutputIterator out)
{
  typename GraphTraits::halfedge_descriptor hb = halfedge(vd,g), done(hb);
  do {
        *out++ = source(hb,g);
        hb = opposite(next(hb,g),g);
  } while(hb!= done);

  return out;
}


template <typename OutputIterator>
OutputIterator
adjacent_vertices_V2(const Polyhedron& g,
                     vertex_descriptor vd,
                     OutputIterator out)
{
  halfedge_around_target_iterator hi, he;

  for(boost::tie(hi, he) = halfedges_around_target(halfedge(vd,g),g); hi != he; ++hi)
  {
    *out++ = source(*hi,g);
  }
  return out;
}


int main(int argc, char** argv)
{
  std::ifstream in((argc>1)?argv[1]:"cube.off");
  Polyhedron P;
  in >> P;
  GraphTraits::vertex_iterator vi = vertices(P).first;
  std::list<vertex_descriptor> V;
  adjacent_vertices_V1(P, *vi, std::back_inserter(V));
  ++vi;
  adjacent_vertices_V2(P, *vi, std::back_inserter(V));
  std::cerr << "done\n";
  return 0;
}
