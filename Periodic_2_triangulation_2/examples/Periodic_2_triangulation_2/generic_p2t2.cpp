#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Generic_P2T2/Periodic_2_Delaunay_triangulation_2_generic.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/Real_timer.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                            K;
typedef K::FT                                                                          FT;
typedef K::Vector_2                                                                    Vector;

typedef typename CGAL::Periodic_2_offset_2                                             Offset;
typedef CGAL::Lattice_2<K>                                                             Lattice;
typedef CGAL::Periodic_2_triangulations_2::internal::Lattice_construct_point_2<K>      CP2;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_base_2<K, Offset, Lattice, CP2> GT;

typedef CGAL::Periodic_2_triangulation_vertex_base_2_generic<GT>                       Vb;
typedef CGAL::Triangulation_face_base_with_info_2<std::pair<bool, CGAL::Color>, GT>    Fbb;
typedef CGAL::Periodic_2_triangulation_face_base_2_generic<GT, Fbb>                    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                                   Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2_generic<GT, Tds>                     GPDT;

typedef GPDT::Point                                                                    Point;

typedef GPDT::Vertex_handle                                                            Vertex_handle;
typedef GPDT::Face_handle                                                              Face_handle;

typedef CGAL::Creator_uniform_2<FT, Point>                                             Creator;
typedef CGAL::Random_points_in_square_2<Point, Creator>                                Point_generator;
typedef CGAL::Random_points_on_circle_2<Point, Creator>                                Point_on_circle_generator;

Lattice generate_random_lattice(CGAL::Random& rnd)
{
  // generate 2 random basis vectors
  std::vector<Point> pts;
  std::copy_n(Point_on_circle_generator(1, rnd), 2, back_inserter(pts));

  std::vector<Vector> dirs;
  for(const Point& pt : pts)
    dirs.emplace_back(CGAL::ORIGIN, pt);

  // random scaling of the directions
  for(Vector& v : dirs)
    v = rnd.get_double(1., 100.) * v;

  std::cout << "Basis vectors:" << std::endl
            << dirs[0] << std::endl
            << dirs[1] << std::endl;

  return Lattice(dirs[0], dirs[1]);
}

std::vector<Point> generate_random_points_in_lattice(const std::size_t n, CGAL::Random& rnd, const Lattice& l)
{
  std::vector<Point> pts;
  pts.reserve(n);

  for(std::size_t i=0; i<n; ++i)
  {
    const double l0 = rnd.get_double();
    const double l1 = rnd.get_double();
    pts.push_back(CGAL::ORIGIN + l0 * l.basis()[0] + l1 * l.basis()[1]);
  }

  return pts;
}

int main(int argc, char** argv)
{
  const std::size_t number_of_points = (argc > 1) ? std::atoi(argv[1]) : 10000;

  const int random_seed = (argc > 2) ? std::atoi(argv[2]) : CGAL::get_default_random().get_int(0, (1 << 30));
  std::cout << "random seed: " << random_seed << std::endl;
  CGAL::Random rnd(random_seed);

  const Lattice l = generate_random_lattice(rnd);
  const std::vector<Point> pts = generate_random_points_in_lattice(number_of_points, rnd, l);

  std::cout << "Inserting " << number_of_points << " random points..." << std::endl;

  CGAL::Real_timer timer;
  timer.start();

  GPDT Tr(l);
  for(const Point& pt : pts)
    Tr.insert(pt);

  std::cout << "Done!" << std::endl;
  std::cout << "Number of vertices: " << Tr.number_of_vertices() << std::endl;
  std::cout << "Number of faces: " << Tr.number_of_faces() << std::endl;
  std::cout << "Is the triangulation 1-cover? " << std::boolalpha << Tr.is_1_cover() << std::endl;

  timer.stop();
  std::cout << "Time: " << timer.time() << " s" << std::endl;

  if(Tr.is_1_cover())
  {
    std::ofstream out("final.off");
    CGAL::write_PD2T2_to_OFF(out, Tr.p2dt2);
  }
  else
  {
    std::ofstream out("final.off");
    CGAL::write_DT2_to_OFF(out, Tr.dt2);
  }

  return EXIT_SUCCESS;
}
