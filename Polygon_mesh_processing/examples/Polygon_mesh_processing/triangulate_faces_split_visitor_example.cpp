#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>

#include <fstream>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                    Point;
typedef CGAL::Surface_mesh<Point>          Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;



class Insert_iterator
{
  typedef boost::unordered_map<face_descriptor,face_descriptor> Container;
  Container& container;
public:

  Insert_iterator(Container &c)
  : container(c) {}

  Insert_iterator&
  operator=(const std::pair<face_descriptor, face_descriptor>& p)
  {
    container[p.second] = p.first;
    return *this;
  }

  Insert_iterator&
  operator*() { return *this; }

  Insert_iterator
  operator++(int) { return *this; }

};


struct Visitor
{
   typedef boost::unordered_map<face_descriptor,face_descriptor> Container;

  Container& container;
  face_descriptor qfd;

  Visitor(Container& container)
    : container(container)
  {}

  void before_subface_creations(face_descriptor fd)
  {
    std::cout << "split : " << fd << " into:" << std::endl;
    Container::iterator it = container.find(fd);
    qfd = it->second;
    container.erase(it);
  }

  void after_subface_created(face_descriptor fd)
  {
    std::cout << "  " << fd;
    container[fd]=qfd;
  }

  void after_subface_creations()
  {
    std::cout << std::endl;
  }
};


int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/P.off";
  std::ifstream input(filename);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  boost::unordered_map<face_descriptor,face_descriptor> t2q;

  Surface_mesh copy;

  copy_face_graph(mesh, copy, CGAL::Emptyset_iterator(), CGAL::Emptyset_iterator(),Insert_iterator(t2q));

  Visitor v(t2q);
  CGAL::Polygon_mesh_processing::triangulate_faces(copy,
                                                   CGAL::Polygon_mesh_processing::parameters::visitor(v));


  for(boost::unordered_map<face_descriptor,face_descriptor>::iterator it = t2q.begin(); it != t2q.end(); ++it){
    std::cout << it->first << "  "  << it->second << std::endl;
  }

  return 0;
}
