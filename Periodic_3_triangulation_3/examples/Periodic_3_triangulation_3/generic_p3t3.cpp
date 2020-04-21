#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#define CGAL_GENERIC_P3T3 // @todo still needed but to remove eventually
#define CGAL_DEBUG_P3T3
#define CGAL_NO_STRUCTURAL_FILTERING

#include <CGAL/internal/Generic_P3T3/Periodic_3_Delaunay_triangulation_3_generic.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/periodic_3_triangulation_3_io.h>
#include <CGAL/IO/Triangulation_off_ostream_3.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                            K;
typedef K::FT                                                                          FT;
typedef K::Vector_3                                                                    Vector;

typedef typename CGAL::Periodic_3_offset_3                                             Offset;
typedef CGAL::Lattice_3<K>                                                             Lattice;
typedef CGAL::Periodic_3_triangulations_3::internal::Lattice_construct_point_3<K>      CP3;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_base_3<K, Offset, Lattice, CP3> GT;

typedef CGAL::Periodic_3_triangulation_vertex_base_3_generic<GT>                       Vb;
typedef CGAL::Periodic_3_triangulation_cell_base_3_generic<GT>                         Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                                   Tds;
typedef CGAL::Periodic_3_Delaunay_triangulation_3_generic<GT, Tds>                     GPDT;

typedef GPDT::Point                                                                    Point;

typedef GPDT::Vertex_handle                                                            Vertex_handle;
typedef GPDT::Cell_handle                                                              Cell_handle;

typedef CGAL::Creator_uniform_3<FT, Point>                                             Creator;
typedef CGAL::Random_points_in_cube_3<Point, Creator>                                  Point_generator;

std::vector<Point> generate_random_points(const std::size_t n, CGAL::Random& rnd)
{
  std::vector<Point> pts;
  pts.reserve(n);

  std::copy_n(Point_generator(3, rnd), n, back_inserter(pts)); // square centered on 0, side 2

  return pts;
}

int main(int argc, char** argv)
{
  const std::size_t number_of_points = (argc > 1) ? std::atoi(argv[1]) : 800;
  std::cout << number_of_points << " random points" << std::endl;
  const int random_seed = (argc > 2) ? std::atoi(argv[2]) : CGAL::get_default_random().get_int(0, (1 << 30));
  std::cout << "random seed: " << random_seed << std::endl;

  CGAL::Random rnd(random_seed);
  const std::vector<Point> pts = generate_random_points(number_of_points, rnd);

  const std::array<Vector, 3> basis = CGAL::make_array(Vector(0, -0.4, 3.95), Vector(-0.3, 0.6, -0.3), Vector(0.4, 0.1, 1.21));
  const CGAL::Lattice_3<K> lattice(basis);

  GPDT Tr(lattice);
  for(const Point& pt : pts)
    Tr.insert(pt);

  std::cout << "1 cover? " << Tr.is_1_cover() << std::endl;
  std::cout << "Number of vertices: " << Tr.number_of_vertices() << std::endl;
//  std::cout << "Number of edges: " << Tr.number_of_edges() << std::endl;
  std::cout << "Number of cells: " << Tr.number_of_cells() << std::endl;

  if(Tr.is_1_cover())
  {
    std::ofstream out("final.off");
    CGAL::write_triangulation_to_off(out, Tr.p3dt3);
  }
  else
  {
    std::ofstream out("final.off");
    CGAL::export_triangulation_3_to_off(out, Tr.dt3);
  }

  return EXIT_SUCCESS;
}
