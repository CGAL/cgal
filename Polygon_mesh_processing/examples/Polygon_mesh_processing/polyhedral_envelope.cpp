#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedral_envelope.h>
#include <CGAL/Surface_mesh.h>

#include <fstream>

int main(int argc, char* argv[])
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point_3;
  typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
  typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;

  typedef CGAL::Polyhedral_envelope<Kernel> Envelope;

  std::ifstream in((argc>1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off"));
  Surface_mesh tmesh;

  in >> tmesh;

  double eps = (argc>2) ? std::stod(std::string(argv[2])) : 0.2;

  Envelope envelope(tmesh, eps);

  int i = (argc>3) ? std::stoi(std::string(argv[3])) : 0;
  int j = (argc>4) ? std::stoi(std::string(argv[4])) : 100;
  int k = (argc>5) ? std::stoi(std::string(argv[5])) : 200;

  if(envelope(tmesh.point(vertex_descriptor(i)),
              tmesh.point(vertex_descriptor(j)),
              tmesh.point(vertex_descriptor(k)))){

    std::cout << "inside polyhedral envelope" << std::endl;
  }

  return 0;
}
