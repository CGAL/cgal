// Benchmark program for the Triangulation_3 package.
//
// Sylvain Pion, 2009.
//
// The data produced by this program is meant to be used
// in the Benchmarks section of the User Manual.
// 
// We measure :
// - construction time
// - point location time function of triangulation size
// - vertex removal time
// - memory usage
//
// And this, for the following 4 configurations :
// - Delaunay
// - Delaunay<Fast_location>
// - Regular
// - Regular<hidden points discarded>
// 
// Notes :
// - points are randomly generated using drand48()
// - weights are zero for regular
//
// TODO (?) :
// - impact of the choice of various kernels
// - impact of the kind of data set ?  More or less degenerate...
// - replace drand48() by CGAL Generators
// - move/document Time_accumulator to CGAL/Profiling_tools (?) 

#define CGAL_TRIANGULATION_3_PROFILING
//#define CGAL_CONCURRENT_TRIANGULATION_3_ADD_TEMPORARY_POINTS_ON_FAR_SPHERE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Memory_sizer.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib> // for drand48

#ifdef SHOW_ITERATIONS
#  undef SHOW_ITERATIONS
#  define SHOW_ITERATIONS "\t( " << iterations << " iterations used)" << endl
#else
#  define SHOW_ITERATIONS endl
#endif

#ifndef BENCH_MIN_TIME
#  define BENCH_MIN_TIME 1 // minimum time for a bench
#endif

using namespace std;
using namespace CGAL;

// Choose the kernel type by defining one of those macros:
// - SC_DOUBLE,
// - EPEC,
// - or EPIC (the default)
#ifdef SC_DOUBLE
   typedef Simple_cartesian<double>                       K;
#elif C_LEDA
#  include <CGAL/leda_rational.h>
#  include <CGAL/Cartesian.h>
   typedef Cartesian<leda_rational> K;
#elif defined(ONLY_STATIC_FILTERS)
   typedef CGAL::internal::Static_filters<CGAL::Simple_cartesian<double> > K;
#elif defined(EPEC)
#  ifdef CGAL_DONT_USE_LAZY_KERNEL
     typedef Epeck K;
#  else // not CGAL_DONT_USE_LAZY_KERNEL
#    ifdef CGAL_USE_LEDA
#      include <CGAL/leda_rational.h>
#      include <CGAL/Cartesian.h>
       typedef Cartesian<leda_rational> SK;
#    else // not CGAL_USE_LEDA
       typedef Simple_cartesian<Gmpq> SK;
#    endif // not CGAL_USE_LEDA
     typedef Lazy_kernel<SK> K;
#  endif // not CGAL_DONT_USE_LAZY_KERNEL
#else // EPIC
  typedef Exact_predicates_inexact_constructions_kernel  K;
#endif // EPIC
typedef Regular_triangulation_euclidean_traits_3<K>    WK;
typedef K::Point_3                                     Point;

#ifdef CGAL_CONCURRENT_TRIANGULATION_3
  typedef CGAL::Spatial_lock_grid_3<
    CGAL::Tag_priority_blocking>                       Lock_ds;

  // Delaunay T3
  typedef CGAL::Triangulation_data_structure_3< 
    CGAL::Triangulation_vertex_base_3<K>, 
    CGAL::Triangulation_cell_base_3<K>, 
    CGAL::Parallel_tag >	                        DT_Tds;
  typedef CGAL::Delaunay_triangulation_3<
    K, DT_Tds, CGAL::Default, Lock_ds>            DT3;
  /*typedef CGAL::Delaunay_triangulation_3<
    K, DT_Tds, Fast_location, Lock_ds>            DT3_FastLoc;*/ // (no parallel fast location for now)

  // Regular T3 with hidden points kept
  typedef CGAL::Triangulation_data_structure_3<
    CGAL::Triangulation_vertex_base_3<WK>, 
    CGAL::Regular_triangulation_cell_base_3<WK>,
    CGAL::Parallel_tag>                           RT_Tds_WithHP;
  typedef CGAL::Regular_triangulation_3<
    WK, RT_Tds_WithHP, Lock_ds>                   RT3_WithHP;

  // Regular T3 with hidden points discarded
  typedef CGAL::Triangulation_data_structure_3<
    CGAL::Triangulation_vertex_base_3<WK>, 
    CGAL::Triangulation_cell_base_3<WK>,
    CGAL::Parallel_tag >	                        RT_Tds_NoHP;
  typedef CGAL::Regular_triangulation_3<
    WK, RT_Tds_NoHP, Lock_ds>                     RT3_NoHP;

#else
  typedef CGAL::Delaunay_triangulation_3<K>                   DT3;
  typedef CGAL::Delaunay_triangulation_3<K, Fast_location>    DT3_FastLoc;
  
  // Regular T3 with hidden points kept
  typedef CGAL::Regular_triangulation_3<WK>                   RT3_WithHP;
  
  // Regular T3 with hidden points discarded
  typedef CGAL::Triangulation_data_structure_3<
    CGAL::Triangulation_vertex_base_3<WK>, 
    CGAL::Triangulation_cell_base_3<WK> >	                    RT_Tds_NoHP;
  typedef CGAL::Regular_triangulation_3<WK, RT_Tds_NoHP>      RT3_NoHP;
#endif

vector<Point> pts, pts2;
Bbox_3 pts_bbox, pts2_bbox;
size_t min_pts = 100;
size_t max_pts = 1000000;

bool input_file_selected = false;
std::ifstream input_file;


class Time_accumulator
{
  double &accumulator;
  Real_timer timer;
public:
  Time_accumulator(double &acc) : accumulator(acc) { timer.reset(); timer.start(); }
  ~Time_accumulator() { timer.stop(); accumulator += timer.time(); }
};

#define drand48 CGAL::get_default_random().get_double

Point rnd_point()
{
  return Point(drand48(), drand48(), drand48());
}

void generate_points()
{
  if (input_file_selected) {
    Point p;
    if (input_file >> p)
    {
      pts.push_back(p);
      pts_bbox = Bbox_3(p.bbox());

      while (input_file >> p)
      {
        pts.push_back(p);
        pts_bbox = pts_bbox + p.bbox();
      }
    }
    cout << " [ Read " << pts.size() << " points from file ] " << endl;
    min_pts = max_pts = pts.size();
  }
  else {
    pts.reserve(max_pts);
    pts2.reserve(max_pts);

    Point p = rnd_point();
    pts.push_back(p);
    pts_bbox = Bbox_3(p.bbox());
    p = rnd_point();
    pts2.push_back(p);
    pts2_bbox = Bbox_3(p.bbox());

    for(size_t i = 1; i < (std::max)(std::size_t(100000), max_pts); ++i) {
      p = rnd_point();
      pts.push_back(p);
      pts_bbox = pts_bbox + p.bbox();
      p = rnd_point();
      pts2.push_back(p);
      pts2_bbox = pts2_bbox + p.bbox();
    }
  }

  cout << "Bounding box = "
    << "[" << pts_bbox.xmin() << ", " << pts_bbox.xmax() << "], "
    << "[" << pts_bbox.ymin() << ", " << pts_bbox.ymax() << "], "
    << "[" << pts_bbox.zmin() << ", " << pts_bbox.zmax() << "]"
    << endl;
}


// Triangulation construction
template < typename Tr >
void benchmark_construction()
{
  cout << "\nTriangulation construction : " << endl;
  cout << "#pts\tTime" << endl;
  size_t mem_size_init = Memory_sizer().virtual_size();
  size_t mem_size = 0;

  for (size_t nb_pts = min_pts; nb_pts <= max_pts; nb_pts *= 10)
  {
    double time = 0;
    size_t iterations = 0;
    do {
      ++iterations;
      Time_accumulator tt(time);
#ifdef CGAL_CONCURRENT_TRIANGULATION_3
      Lock_ds locking_ds(pts_bbox, 50);
      Tr tr(pts.begin(), pts.begin() + nb_pts, &locking_ds);
#else
      Tr tr(pts.begin(), pts.begin() + nb_pts);
#endif
      mem_size = Memory_sizer().virtual_size();
      // cout << "#vertices = " << tr.number_of_vertices() << endl;
      // cout << "#cells = " << tr.number_of_cells() << endl;
    } while (time < BENCH_MIN_TIME);
    cout << nb_pts << "\t" << time/iterations << SHOW_ITERATIONS;
  }
  cout << "\nMemory usage : " << (mem_size - mem_size_init)*1./max_pts << " Bytes/Point"
       << " (observed for the largest data set, and may be unreliable)" << endl;
}


// Point location
template < typename Tr >
void benchmark_location()
{
  cout << "\nPoint location : " << endl;
  cout << "#pts\tTime" << endl;
  for (size_t nb_pts = min_pts; nb_pts <= max_pts; nb_pts *= 10)
  {
#ifdef CGAL_CONCURRENT_TRIANGULATION_3
    Lock_ds locking_ds(pts_bbox, 50);
    Tr tr(pts.begin(), pts.begin() + nb_pts, &locking_ds);
#else
    Tr tr(pts.begin(), pts.begin() + nb_pts);
#endif
    double time = 0;
    size_t iterations = 0;
    do {
      ++iterations;
      Time_accumulator tt(time);
      // We do chunks of    100000 point locations at once.
      for(size_t i = 0; i < 100000; ++i)
        tr.locate(pts2[i]);
    } while (time < BENCH_MIN_TIME);
    cout << nb_pts << "\t" << (time/iterations)/100000 << SHOW_ITERATIONS;
  }
}


// Vertex removal
template < typename Tr >
void benchmark_remove()
{
  typedef typename Tr::Vertex_handle     Vertex_handle;
  typedef typename Tr::Vertex_iterator   Vertex_iterator;

  cout << "\nVertex removal : " << endl;
  cout << "#pts\tTime\tTime/removal" << endl;
  size_t nb_pts = 1000000; // only one size of triangulation hard-coded.
  const size_t NUM_VERTICES_TO_REMOVE = 100000;
  double time = 0;
  size_t iterations = 0;

  if (nb_pts > max_pts)
  {
    std::cerr << "ERROR: nb_pts > max_pts. Cancelling..." << std::endl;
    return;
  }

  do {
#ifdef CGAL_CONCURRENT_TRIANGULATION_3
    Lock_ds locking_ds(pts_bbox, 50);
    Tr tr(pts.begin(), pts.begin() + nb_pts, &locking_ds);
#else
    Tr tr(pts.begin(), pts.begin() + nb_pts);
#endif
    vector<Vertex_handle> vhs;
    for (Vertex_iterator vit = tr.finite_vertices_begin(), end = tr.finite_vertices_end();
         vit != end; ++vit)
      vhs.push_back(vit);

    Time_accumulator tt(time);
    tr.remove(&vhs[0], &vhs[NUM_VERTICES_TO_REMOVE - 1]);
    ++iterations;
				//std::cout<<"\b\b\b\b\b\b"<<i<<std::flush;
				//tr.is_valid();
			}
  } while (time < BENCH_MIN_TIME);

  cout << NUM_VERTICES_TO_REMOVE << "\t" 
        << (time/iterations) << "\t" 
        << (time/iterations)/NUM_VERTICES_TO_REMOVE << SHOW_ITERATIONS;
}


template < typename Tr >
void do_benchmarks(string name)
{
  cout << "\n\nBenchmarking configuration : " << name << endl;
  // tbb::task_scheduler_init tbb_init(10); // Set number of threads
  benchmark_construction<Tr>();
  if (input_file_selected)
    return;
  benchmark_location<Tr>();
  benchmark_remove<Tr>();
}

int main(int argc, char **argv)
{
  if (argc >= 2) {
    input_file.open(argv[1], std::ios::in);
    if (input_file.is_open())
      input_file_selected = true;
    else {
      input_file_selected = false;
      max_pts = atoi(argv[1]);
    }
  }

  cout << "Usage : " << argv[0] << " [filename]"
       << " [max_points = " << max_pts << ", and please use a power of 10]" << endl;
  cout << "Benchmarking the Triangulation_3 package for ";
  if (input_file_selected)
    cout << "data file : " << argv[1] << endl;
  else
    cout << "up to " << max_pts << " random points." << endl;

  cout.precision(3);

  generate_points();

  cout << "\nProcessor : "
       << ((sizeof(void*)==4) ? 32 : (sizeof(void*)==8) ? 64 : -1) << " bits\n";
  cout << "Kernel typeid: " << typeid(K).name() << "\n";

  do_benchmarks<DT3>("Delaunay  [Compact_location]");
  if (input_file_selected)
    return 0;
#ifndef CGAL_CONCURRENT_TRIANGULATION_3
  do_benchmarks<DT3_FastLoc>("Delaunay with Fast_location");
#endif
  do_benchmarks<RT3_WithHP>("Regular  [with hidden points kept, except there's none in the data sets]");
  do_benchmarks<RT3_NoHP>("Regular with hidden points discarded");
}
