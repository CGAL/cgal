#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#define CGAL_GENERIC_P2T2 // @todo still needed but to remove eventually

// #define CGAL_DEBUG_P2T2

#include <CGAL/internal/Generic_P2T2/Periodic_2_Delaunay_triangulation_2_generic.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                              K;
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

typedef CGAL::Creator_uniform_2<FT, Point>                        Creator;
typedef CGAL::Random_points_in_square_2<Point, Creator>           Point_generator;
typedef CGAL::Real_timer         Timer;


std::vector<Point> generate_random_points(std::size_t n, CGAL::Random& rnd)
{
  std::vector<Point> pts;
  pts.reserve(n);

  std::copy_n(Point_generator(1, rnd), n, back_inserter(pts)); // square centered on 0, side 2

  return pts;
}

int main(int argc, char** argv)
{
  std::ofstream log = std::ofstream("log.txt");
  const std::size_t number_of_points = (argc > 1) ? std::atoi(argv[1]) : 800;
  log << number_of_points << " random points" << std::endl;
  const std::size_t ntrials = (argc > 2) ? std::atoi(argv[2]) : 100;
  log << ntrials << " trials" << std::endl;
  std::array<Vector, 2> basis;
  if (argc > 6) {
    basis = CGAL::make_array(Vector(std::atof(argv[3]), std::atof(argv[4])), Vector(std::atof(argv[5]), std::atof(argv[6])));
  } else {
    basis = CGAL::make_array(Vector(1, 0), Vector(0, 1));
  }
  log << "basis: " << basis[0] << " ; " << basis[1] << std::endl;
  const int random_seed = (argc > 7) ? std::atoi(argv[7]) : CGAL::get_default_random().get_int(0, (1 << 30));
  log << "random seed: " << random_seed << std::endl;
  const CGAL::Lattice_2<K> lattice0(basis);
  log << "initialization phase1 conversion phase2 total" << std::endl;

  Timer timer = Timer();
  double t_lattice=0, t_phase1=0, t_switch=0, t_phase2=0, t_total=0;
  for (int trial = 0; trial < ntrials; ++trial) {  
    CGAL::Random rnd(random_seed + trial);
    std::vector<Point> pts = generate_random_points(number_of_points, rnd);
    std::transform(pts.begin(), pts.end(), pts.begin(),
                   [lattice0](const Point& pt) -> Point { return lattice0.construct_canonical_point(pt); });

    double t, ts, t2, t3;
    timer.reset();
    timer.start();
    const CGAL::Lattice_2<K> lattice(basis);
    GPDT Tr(lattice);
    t = timer.time();
    t_lattice += t;
    log << t << " ";

    // Tr.insert(pts.begin(), pts.end());
    bool switch_occurred = false;
    for(const Point& pt : pts)
    {
      // if(lattice.is_in_scaled_domain(pt))
      ts = timer.time();
      Tr.insert(pt);
      if (!switch_occurred && Tr.is_1_cover()) {
        t2 = timer.time();
        t_phase1 += ts - t;
        t_switch += t2 - ts;
        log << ts - t << " ";
        log << t2 - ts << " ";
        switch_occurred = true;
      }
    }
    t3 = timer.time();
    timer.stop();
    if (Tr.is_1_cover()) {
      t_phase2 += t3-t2;
    } else {
      t_phase1 += t3 - t;
    }
    t_total += t3;
    log << t3-t2 << " " << t3 << std::endl;
  }
  log << "Average: " << t_lattice/ntrials << " " << t_phase1/ntrials << " " << t_switch/ntrials << " " << t_phase2/ntrials << " " << t_total/ntrials << std::endl;

  return EXIT_SUCCESS;
}
