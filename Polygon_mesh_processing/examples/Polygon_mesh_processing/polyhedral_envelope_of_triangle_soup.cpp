#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedral_envelope.h>
#include <CGAL/IO/OFF.h>

#include <vector>
#include <fstream>

int main(int argc, char* argv[])
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point_3;

  typedef CGAL::Polyhedral_envelope<Kernel> Envelope;

  std::ifstream in((argc>1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off"));
  double eps = (argc>2) ? std::stod(std::string(argv[2])) : 0.2;


  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > polygons;

  CGAL::IO::read_OFF(in, points, polygons);

  Envelope envelope(points, polygons, eps);

  int i = (argc>3) ? std::stoi(std::string(argv[3])) : 0;
  int j = (argc>4) ? std::stoi(std::string(argv[4])) : 100;
  int k = (argc>5) ? std::stoi(std::string(argv[5])) : 200;

  if (envelope(points[i], points[j],points[k]))
  {
    std::cout << "inside polyhedral envelope" << std::endl;
  }

  return 0;
}
