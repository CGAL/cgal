#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#define CGAL_GENERIC_P3T3 // @todo still needed but to remove eventually
// #define CGAL_DEBUG_P3T3
// #define CGAL_NO_STRUCTURAL_FILTERING

#include <CGAL/internal/Generic_P3T3/Periodic_3_Delaunay_triangulation_3_generic.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/periodic_3_triangulation_3_io.h>
#include <CGAL/IO/Triangulation_off_ostream_3.h>

#include <CGAL/Real_timer.h>

#include <iostream>
#include <fstream>

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
typedef CGAL::Random_points_in_cube_3<Point, Creator>                                  Point_in_cube_generator;
typedef CGAL::Random_points_on_sphere_3<Point, Creator>                                Point_on_sphere_generator;

// Lattice and Volume / cubed systole
Lattice generate_random_lattice(CGAL::Random& rnd)
{
  // generate 3 random basis vectors
  std::vector<Point> pts;
  std::copy_n(Point_on_sphere_generator(1, rnd), 3, back_inserter(pts));

  std::vector<Vector> dirs;
  for(const Point& pt : pts)
    dirs.emplace_back(CGAL::ORIGIN, pt);

  // random scaling of the directions
  for(Vector& v : dirs)
    v = rnd.get_double(1., 100.) * v;

  std::cout << "Basis vectors:" << std::endl
            << dirs[0] << std::endl
            << dirs[1] << std::endl
            << dirs[2] << std::endl;

  return Lattice(dirs[0], dirs[1], dirs[2]);
}

std::vector<Point> generate_random_points(const std::size_t n, CGAL::Random& rnd)
{
  std::vector<Point> pts;
  pts.reserve(n);

  std::copy_n(Point_in_cube_generator(100, rnd), n, back_inserter(pts)); // square centered on 0, side 200

  return pts;
}

std::vector<Point> generate_random_points_in_lattice(const std::size_t n, CGAL::Random& rnd, const Lattice& l)
{
  std::vector<Point> pts;
  pts.reserve(n);

  for(std::size_t i=0; i<n; ++i)
  {
    const double l0 = rnd.get_double();
    const double l1 = rnd.get_double();
    const double l2 = rnd.get_double();
    pts.push_back(CGAL::ORIGIN + l0 * l.basis()[0] + l1 * l.basis()[1] + l2 * l.basis()[2]);
  }

  return pts;
}


int main(int argc, char** argv)
{
  const std::size_t number_of_points = (argc > 1) ? std::atoi(argv[1]) : 1000000;
//  std::cout << number_of_points << " random points" << std::endl;

  const int random_seed = (argc > 2) ? std::atoi(argv[2]) : CGAL::get_default_random().get_int(0, (1 << 30));
//  std::cout << "random seed: " << random_seed << std::endl;
  CGAL::Random rnd(random_seed);

  std::ofstream out("/home/mrouxell/log_200_runs_1000.txt");
  out.precision(17);

  for(;;)
  {
    Lattice l = generate_random_lattice(rnd);

    const FT vol = CGAL::determinant(l.basis()[0][0], l.basis()[0][1], l.basis()[0][2],
                                     l.basis()[1][0], l.basis()[1][1], l.basis()[1][2],
                                     l.basis()[2][0], l.basis()[2][1], l.basis()[2][2]);
    const FT cubed_sys = l.systole_sq_length() * l.systole_sq_length() * l.systole_sq_length();

    int number_of_runs = 200;
//    std::cout << "number_of_runs: " << number_of_runs << std::endl;

    std::cout << "vol/sys3: " << vol/cubed_sys << std::endl;

    bool all_good = true;

    int switch_nv = 0;
    double time = 0.;
    for(int i=0; i<number_of_runs; ++i)
    {
      const std::vector<Point> pts = generate_random_points_in_lattice(number_of_points, rnd, l);

      CGAL::Real_timer timer;
      timer.start();

      GPDT Tr(l);
      for(const Point& pt : pts)
      {
        Tr.insert(pt);
        if(Tr.is_1_cover())
          break;
      }

//      std::cout << "Number of vertices: " << Tr.number_of_vertices() << std::endl;
//      std::cout << "Number of cells: " << Tr.number_of_cells() << std::endl;
      std::cout << "1 cover? " << Tr.is_1_cover() << std::endl;

      if(!Tr.is_1_cover())
      {
        all_good = false;
        break; // break or --number_of_runs?
      }

      timer.stop();
      time += timer.time();
      switch_nv += Tr.switch_nv;

      std::cout << "Switch NV: " << Tr.switch_nv << std::endl;
      std::cout << "Time: " << timer.time() << std::endl;
    }

    std::cout << "Computed in: " << time << " s" << std::endl;

    if(!all_good)
      continue;

    switch_nv /= number_of_runs;
    time /= number_of_runs;

    // only print if we managed to convert to 1 cover within the allowed amount of vertices (i.e. time)

    out << l.basis()[0] << " ";
    out << l.basis()[1] << " ";
    out << l.basis()[2] << " ";

    out << CGAL::abs(vol / cubed_sys) << " ";
    out << switch_nv << " ";

    out << time << std::endl;

    //    if(Tr.is_1_cover())
    //    {
    //      std::ofstream out("final.off");
//      CGAL::write_triangulation_to_off(out, Tr.p3dt3);
//    }
//    else
//    {
//      std::ofstream out("final.off");
//      CGAL::export_triangulation_3_to_off(out, Tr.dt3);
//    }
  }

  return EXIT_SUCCESS;
}
