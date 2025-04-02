// #define DOUBLE_2D_SNAP_VERBOSE

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Float_snap_rounding_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Curve_2                                       Segment_2;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::Vector_2                                        Vector_2;
typedef std::vector<Point_2 >                                   Polyline_2;
typedef Kernel::FT                                              FT;

Point_2 random_point(CGAL::Random& r, double m = -1, double M = 1)
{
  return Point_2(FT(r.get_double(m, M)), FT(r.get_double(m, M)));
}

Point_2 random_noise_point(CGAL::Random& r, double m, double M, Vector_2 v)
{
  return Point_2(FT(r.get_double(m, M) + CGAL::to_double(v.x())), FT(r.get_double(m, M)  + CGAL::to_double(v.y())));
}

void test_fully_random(CGAL::Random &r, size_t nb_segments){
  std::cout << "Test fully random" << std::endl;
  std::vector<Segment_2> segs;
  std::vector<Segment_2> out;

  for(size_t i=0; i!=nb_segments; ++i)
    segs.emplace_back(random_point(r), random_point(r));

  CGAL::Real_timer t;
  t.start();
  CGAL::compute_snap_subcurves_2(segs.begin(), segs.end(), out);
  t.stop();

  std::cout << t.time() << "sec" << std::endl;
  assert(!CGAL::do_curves_intersect(out.begin(), out.end()));
}

void test_almost_indentical_segments(CGAL::Random &r, size_t nb_segments, Vector_2 source, Vector_2 target){
  // Difficult test
  std::cout << "Test with almost identical segments from " << source << " to " << target << std::endl;
  std::vector<Segment_2> segs;
  std::vector<Segment_2> out;

  double eps=std::pow(2,-45);
  for(size_t i=0; i!=nb_segments; ++i){
    segs.emplace_back(random_noise_point(r, -eps, eps, source), random_noise_point(r, -eps, eps, target));
  }

  CGAL::Real_timer t;
  t.start();
  CGAL::compute_snap_subcurves_2(segs.begin(), segs.end(), out);
  t.stop();

  std::cout << t.time() << "sec" << std::endl;
  assert(!CGAL::do_curves_intersect(out.begin(), out.end()));
}

void test_multi_almost_indentical_segments(CGAL::Random &r, size_t nb_segments){
  for(double x1=-1; x1<=1; ++x1)
    for(double y1=-1; y1<=1; ++y1)
      for(double x2=x1; x2<=1; ++x2)
        for(double y2=y1; y2<=1; ++y2)
          if(Vector_2(x1,y1)!=Vector_2(x2,y2))
            test_almost_indentical_segments(r, nb_segments, Vector_2(x1,y1), Vector_2(x2,y2));
}

int main(int argc,char *argv[])
{
  CGAL::Random rp;
  CGAL::Random r(argc==1?rp.get_seed():std::stoi(argv[1]));
  std::cout << "random seed = " << r.get_seed() << std::endl;
  std::cout << std::setprecision(17);
  test_data();
  test_fully_random(r,2000);
  test_multi_almost_indentical_segments(r,200);
  return(0);
}
