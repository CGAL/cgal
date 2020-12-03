#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedral_envelope.h>
#include <CGAL/Surface_mesh.h>

#include <fstream>

int main(int argc, char* argv[])
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point_3;
  typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
  typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::Polyhedral_envelope<Kernel> Envelope;

  std::ifstream in((argc>1) ? argv[1] : "data/patch.off");
  std::ifstream in2((argc>2) ? argv[2] : "data/tentative.off");
  Surface_mesh tmesh, test;

  in >> tmesh;
  in2 >> test;

  double eps = (argc>3) ? std::stod(std::string(argv[3])) : 0.2;

  Envelope envelope;
  envelope = Envelope(tmesh, eps);

  for(auto fd : faces(test)){
    halfedge_descriptor hd = halfedge(fd, test);
    Point_3 p = test.point(target(hd,test));
    Point_3 q = test.point(target(next(hd,test),test));
    Point_3 r = test.point(source(hd,test));
    if(envelope(p,q,r)){
      std::cout << "inside polyhedral envelope" << std::endl;
    }
  }
  std::cout << "done" << std::endl;
  return 0;
}
