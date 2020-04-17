#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/tetrahedral_remeshing_io.h>

#include <CGAL/Random.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;
typedef Remeshing_triangulation::Point          Point;
typedef Remeshing_triangulation::Cell_handle    Cell_handle;


int main(int argc, char* argv[])
{
  const unsigned int nbv = (argc > 1) ? atoi(argv[1]) : 100;

  CGAL::Random rng;
  std::cout << "CGAL Random seed = " << CGAL::get_default_random().get_seed() << std::endl;

  std::vector<Point> points;
  while (points.size() < nbv)
  {
    Point p(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.));
    points.push_back(p);
  }


  Remeshing_triangulation tr(points.begin(), points.end());

  for (Cell_handle c : tr.finite_cell_handles())
    c->set_subdomain_index(1);

  std::cout << "save_ascii_triangulation : ";
  std::cout.flush();
  std::ofstream out1("remeshing_triangulation.ascii.cgal",
                     std::ios_base::out);
  bool ok = CGAL::save_ascii_triangulation(out1, tr);
  out1.close();
  assert(ok);
  std::cout << "done." << std::endl;

  Remeshing_triangulation tr1;
  std::cout << "load_triangulation (ascii) : ";
  std::cout.flush();
  std::ifstream in1("remeshing_triangulation.ascii.cgal",
                    std::ios_base::in);
  ok = CGAL::load_triangulation(in1, tr1);
  assert(ok);
  std::cout << "done." << std::endl;

  std::cout << "save_binary_triangulation : ";
  std::cout.flush();
  std::ofstream out2("remeshing_triangulation.binary.cgal",
                    std::ios_base::out | std::ios_base::binary);
  ok = CGAL::save_binary_triangulation(out2, tr);
  out2.close();
  assert(ok);
  std::cout << "done." << std::endl;

  Remeshing_triangulation tr2;
  std::cout << "load_triangulation (binary) : ";
  std::cout.flush();
  std::ifstream in2("remeshing_triangulation.binary.cgal",
                    std::ios_base::in | std::ios_base::binary);
  ok = CGAL::load_triangulation(in2, tr2);
  assert(ok);
  std::cout << "done." << std::endl;

  return EXIT_SUCCESS;
}
