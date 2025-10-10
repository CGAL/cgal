#define BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Float_snap_rounding_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

#include <CGAL/Boolean_set_operations_2.h>

#include <CGAL/Cartesian_converter.h>

#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef CGAL::Float_snap_rounding_traits_2<Kernel>              Traits_2;
typedef Kernel::Segment_2                                       Segment_2;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::Vector_2                                        Vector_2;
typedef std::vector<Point_2 >                                   Polyline_2;
typedef std::vector<Polyline_2>                                 Polyline_range_2;
typedef Kernel::FT                                              FT;

typedef CGAL::Arr_segment_traits_2<Kernel>                      Arr_segment_traits_2;
typedef Arr_segment_traits_2::Curve_2                           Curve_2;

typedef CGAL::Snap_rounding_traits_2<Kernel>     SnapTraits;

typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef std::vector<Polygon_with_holes_2>                   Pwh_vec_2;

typedef CGAL::Cartesian<double>                                        Naive;
typedef CGAL::Cartesian_converter<Kernel,Naive>                         EK_to_IK;
typedef CGAL::Cartesian_converter<Naive, Kernel>                         IK_to_EK;

//Biggest double with ulp smaller than an integer
constexpr double maxFloat = std::pow(2,23);
constexpr double maxDouble = std::pow(2,52);

Segment_2 scale_segment(const Segment_2 &s){
  auto scale_value=[](FT x){
    return maxDouble * CGAL::to_double(x);
  };
  return Segment_2(Point_2(scale_value(s.start().x()), scale_value(s.start().y())),
                   Point_2(scale_value(s.end().x()), scale_value(s.end().y())));
}

Point_2 random_point(CGAL::Random& r, double m = -1, double M = 1)
{
  return Point_2(FT(r.get_double(m, M)), FT(r.get_double(m, M)));
}

Point_2 random_noise_point(CGAL::Random& r, double m, double M, Vector_2 v)
{
  return Point_2(FT(r.get_double(m, M) + CGAL::to_double(v.x())), FT(r.get_double(m, M)  + CGAL::to_double(v.y())));
}

template<class Iterator, class OutRange>
void compute_subcurves_and_naive_snap(Iterator begin, Iterator end, OutRange &out){
  EK_to_IK to_inexact;
  IK_to_EK to_exact;
  std::vector<Curve_2> segs;
  for(Iterator it=begin; it!=end; ++it)
    segs.emplace_back(*it);
  CGAL::compute_subcurves(segs.begin(), segs.end(), std::back_inserter(out));
  //Naive snap
  OutRange out2;
  for(Segment_2 &seg: out){
    Segment_2 s = to_exact(to_inexact(seg));
    if(!s.is_degenerate())
      out2.emplace_back(s);
  }
  std::swap(out,out2);
}

void test(const std::vector<Segment_2> &segs){
  std::vector<Segment_2> out;
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
  CGAL::Real_timer t;
  std::cout << "Input size: " << segs.size() << std::endl;
  t.start();
  compute_subcurves_and_naive_snap(segs.begin(), segs.end(), out);
  t.stop();
  std::cout << "Naive snap size: " << out.size() << " ,running time: " << t.time() << " and output do intersect: " << CGAL::do_curves_intersect(out.begin(), out.end()) <<std::endl;
  out.clear();
  t.reset();
  t.start();
#endif
  CGAL::compute_snapped_subcurves_2(segs.begin(), segs.end(), out);
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
  t.stop();
  std::cout << "Formal snap size: " << out.size() << " ,running time: " << t.time() << "\n" << std::endl;
#endif
  assert(!CGAL::do_curves_intersect(out.begin(), out.end()));
}

void test_fully_random(CGAL::Random &r, size_t nb_segments){
  std::cout << "Test fully random" << std::endl;
  std::vector<Segment_2> segs;
  for(size_t i=0; i!=nb_segments; ++i)
    segs.emplace_back(random_point(r), random_point(r));

  test(segs);
  // CGAL::Real_timer t;
  // Polyline_range_2 output_list;
  // t.start();
  // CGAL::snap_rounding_2<SnapTraits>(segs.begin(), segs.end(), output_list, 1./maxFloat);
  // t.stop();
  // std::cout << "Running time with integer 2D snap: " << t.time() << "sec" << std::endl;
  // t.reset();
}

void test_almost_indentical_segments(CGAL::Random &r, size_t nb_segments, Vector_2 source, Vector_2 target){
  // Difficult test
  std::cout << "Test with almost identical segments from " << source << " to " << target << std::endl;
  std::vector<Segment_2> segs;
  double eps=std::pow(2,-40);
  for(size_t i=0; i!=nb_segments; ++i){
    segs.emplace_back(random_noise_point(r, -eps, eps, source), random_noise_point(r, -eps, eps, target));
  }

  test(segs);
}

void test_iterative_square_intersection(CGAL::Random &r, size_t nb_iterations){
  auto add_random_rotated_square=[&](std::vector<Segment_2> &segs){
    double theta=r.get_double(0, CGAL_PI/2);
    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);
    Point_2 a( cos_t, sin_t);
    Point_2 b( sin_t,-cos_t);
    Point_2 c(-cos_t,-sin_t);
    Point_2 d(-sin_t, cos_t);
    segs.emplace_back(a,b);
    segs.emplace_back(b,c);
    segs.emplace_back(c,d);
    segs.emplace_back(d,a);
  };

  std::vector<Segment_2> segs;
  std::vector<Segment_2> out;
  std::vector<Curve_2> arr_segs;

  CGAL::Real_timer t;
  for(size_t i=0; i<nb_iterations; ++i){
    std::cout << "Iterations " << i << std::endl;
    out.clear();
    arr_segs.clear();

    for(int j=0; j<5; ++j)
      add_random_rotated_square(segs);

    test(segs);
    CGAL::compute_snapped_subcurves_2(segs.begin(), segs.end(), out);
    assert(!CGAL::do_curves_intersect(out.begin(), out.end()));

    segs.clear();
    segs.insert(segs.begin(), out.begin(), out.end());
  }


}

void test_iterative_square_intersection(CGAL::Random &r, size_t nb_iterations){
  auto random_rotated_square=[&](){
    double theta=r.get_double(0, CGAL_PI/2);
    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);
    Polygon_2 P;
    P.push_back( cos_t, sin_t);
    P.push_back( sin_t,-cos_t);
    P.push_back(-cos_t,-sin_t);
    P.push_back(-sin_t, cos_t);
    return P;
  };

  Pwh_vec_2 scene, out_intersection;

  CGAL::intersection(random_rotated_square, random_rotated_square, scene);

  for(size_t i=0; i<nb_iterations; ++i){
    CGAL::intersection(random_rotated_square, random_rotated_square, out_intersection);

  }


}

void test_multi_almost_indentical_segments(CGAL::Random &r, size_t nb_segments){
  for(double x1=-1; x1<=1; ++x1)
    for(double y1=-1; y1<=1; ++y1)
      for(double x2=x1; x2<=1; ++x2)
        for(double y2=y1; y2<=1; ++y2)
          if(Vector_2(x1,y1)!=Vector_2(x2,y2))
            test_almost_indentical_segments(r, nb_segments, Vector_2(x1,y1), Vector_2(x2,y2));
}

void fix_test(){
  std::cout << "Fix tests" << std::endl;

  std::vector< Segment_2 > segs;
  FT e(std::pow(2, -60));
  segs.emplace_back(Point_2(0, 1+e), Point_2(2, 1));
  segs.emplace_back(Point_2(1, 1), Point_2(1, 0));

  test(segs);
  segs.clear();

  segs.emplace_back(Point_2(1-e, 1), Point_2(-1-e, -1+2*e));
  segs.emplace_back(Point_2(e/2, e/2), Point_2(1, -1));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(2, 0));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(-2+e, -4));
  segs.emplace_back(Point_2(-2, 2), Point_2(2, 2));
  segs.emplace_back(Point_2(7, 7), Point_2(7+e, 7+e));
  segs.emplace_back(Point_2(5, 7-e), Point_2(9, 7-e));
  test(segs);
}

int main(int argc,char *argv[])
{
  CGAL::Random rp;
  CGAL::Random r(argc==1?rp.get_seed():std::stoi(argv[1]));
  std::cout << "random seed = " << r.get_seed() << std::endl;
  std::cout << std::setprecision(17);
  fix_test();
  test_fully_random(r,1000);
  test_multi_almost_indentical_segments(r,100);
  test_iterative_square_intersection(r,500);
  return(0);
}
