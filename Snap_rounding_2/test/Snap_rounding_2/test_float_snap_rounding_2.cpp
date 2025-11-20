#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Float_snap_rounding_traits_2.h>
#include <CGAL/Float_snap_rounding_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

#include <CGAL/Boolean_set_operations_2.h>

#include <CGAL/Random.h>

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

typedef CGAL::Polygon_2<Kernel>                                 Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                      Polygon_with_holes_2;
typedef std::vector<Polygon_with_holes_2>                       Pwh_vec_2;

#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
#include <CGAL/Real_timer.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<double>                          Naive;
typedef CGAL::Cartesian_converter<Kernel,Naive>                 EK_to_IK;
typedef CGAL::Cartesian_converter<Naive, Kernel>                IK_to_EK;
#endif

#ifdef COMPARE_WITH_INTEGER_SNAP_ROUNDING_2
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>
typedef CGAL::Snap_rounding_traits_2<Kernel>                     SnapTraits;
#endif

//Biggest double with ulp smaller than an integer
const double maxFloat = std::pow(2,23);
const double maxDouble = std::pow(2,52);

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

#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
template<class Iterator, class OutRange>
void compute_subcurves_and_naive_snap(Iterator begin, Iterator end, OutRange &out){
  EK_to_IK to_inexact;
  IK_to_EK to_exact;
  std::vector<Curve_2> segs;
  for(Iterator it=begin; it!=end; ++it)
    if(it->source()!=it->target())
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
#endif

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
  CGAL::compute_snapped_subcurves_2(segs.begin(), segs.end(), std::back_inserter(out));
  std::vector<Segment_2> out2(out.begin(), out.end());
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
  t.stop();
  std::cout << "Formal snap size: " << out.size() << " ,running time: " << t.time() << std::endl;
#endif
  assert(!CGAL::do_curves_intersect(out.begin(), out.end()));
#ifdef COMPARE_WITH_INTEGER_SNAP_ROUNDING_2
  Polyline_range_2 output_list;
  t.reset();
  t.start();
  CGAL::snap_rounding_2<SnapTraits>(segs.begin(), segs.end(), std::back_inserter(output_list), 1./maxDouble);
  t.stop();
  std::cout << "Running time with integer 2D snap (scaled at 10^15): " << t.time() << std::endl;
#endif
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
  std::cout << "\n";
#endif
}

void test_fully_random(CGAL::Random &r, size_t nb_segments){
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
  std::cout << "Test fully random" << std::endl;
#endif
  std::vector<Segment_2> segs;
  for(size_t i=0; i!=nb_segments; ++i)
    segs.emplace_back(random_point(r), random_point(r));

  test(segs);
}

void test_random_polygons(CGAL::Random &r, size_t nb_polygons, size_t nb_pts){
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
  std::cout << "Test random polygons" << std::endl;
#endif
  std::vector<Polygon_2> polygons;
  for(size_t i=0; i!=nb_polygons; ++i){
    Polygon_2 poly;
    for(size_t j=0; j!=nb_pts; ++j)
      poly.push_back(Point_2(random_point(r), random_point(r)));
    polygons.emplace_back(poly);
  }

  std::vector<Polygon_2> out;
  CGAL::compute_snapped_polygons_2(polygons.begin(), polygons.end(), std::back_inserter(out));

  std::vector<Segment_2> segs;
  for(const Polygon_2 &poly: out)
    for(size_t i=1; i<poly.size(); ++i)
      segs.emplace_back(poly[i-1],poly[i]);
  assert(!CGAL::do_curves_intersect(segs.begin(), segs.end()));
}

void test_almost_indentical_segments(CGAL::Random &r, size_t nb_segments, Vector_2 source, Vector_2 target){
  // Difficult test
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
  std::cout << "Test with almost identical segments from " << source << " to " << target << std::endl;
#endif
  std::vector<Segment_2> segs;
  double eps=std::pow(2,-40);
  for(size_t i=0; i!=nb_segments; ++i){
    segs.emplace_back(random_noise_point(r, -eps, eps, source), random_noise_point(r, -eps, eps, target));
  }

  test(segs);
}

void test_iterative_square_intersection(CGAL::Random &r, size_t nb_iterations){
  auto random_rotated_square=[&](){
    double theta=r.get_double(0, CGAL_PI/2);
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
    std::cout << "Angle: " << theta << std::endl;
#endif
    FT cos_t(std::cos(theta));
    FT sin_t(std::sin(theta));
    Polygon_2 P;
    P.push_back(Point_2( cos_t, sin_t));
    P.push_back(Point_2(-sin_t, cos_t));
    P.push_back(Point_2(-cos_t,-sin_t));
    P.push_back(Point_2( sin_t,-cos_t));
    return P;
  };

  Polygon_2 scene=random_rotated_square();
  Polygon_2 snap_scene;
  Pwh_vec_2 out_intersection;

  for(size_t i=0; i<nb_iterations; ++i){
    out_intersection.clear();
    CGAL::intersection(random_rotated_square(), scene, std::back_inserter(out_intersection));
    assert(out_intersection.size()==1 && out_intersection[0].number_of_holes()==0);
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
    CGAL::Real_timer t;
    t.start();
#endif
    compute_snapped_polygon_2(out_intersection[0].outer_boundary(), snap_scene);
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
    t.stop();
    std::cout << "Iteration " << i << std::endl;
    std::cout << "Polygon size: " << out_intersection[0].outer_boundary().size()
              << " , Snapped polygon size: " << snap_scene.size() << std::endl;
    std::cout << "is convex: " << snap_scene.is_convex() << std::endl;
    std::cout << "Running time: " << t.time() << std::endl;
#endif
    scene=snap_scene;
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
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
  std::cout << "Fix tests" << std::endl;
#endif
  std::vector< Segment_2 > segs;
  FT e(std::pow(2, -60));

  // Disjoint segments
  segs.emplace_back(Point_2(0, 1), Point_2(1, 1));
  segs.emplace_back(Point_2(1, 2), Point_2(2, 2));
  test(segs);
  segs.clear();

  // Degenerate segment
  segs.emplace_back(Point_2(0, 1), Point_2(2, 1));
  segs.emplace_back(Point_2(1, 0), Point_2(1, 0));
  test(segs);
  segs.clear();

  // Collinear segments
  segs.emplace_back(Point_2(0, 0), Point_2(4, 0));
  segs.emplace_back(Point_2(2, 0), Point_2(6, 0));
  test(segs);
  segs.clear();

  // Duplicate segment
  segs.emplace_back(Point_2(0, 0), Point_2(4, 0));
  segs.emplace_back(Point_2(0, 0), Point_2(4, 0));
  test(segs);
  segs.clear();

  // Polyline
  segs.emplace_back(Point_2(0, 0), Point_2(2, 0));
  segs.emplace_back(Point_2(2, 0), Point_2(1, 1));
  segs.emplace_back(Point_2(1, 1), Point_2(3, 1));
  test(segs);
  segs.clear();

  // Horizontal intersection when round
  segs.emplace_back(Point_2(0, 1+e), Point_2(2, 1));
  segs.emplace_back(Point_2(1, 1), Point_2(1, 0));
  test(segs);
  segs.clear();

  // Vertical intersection when round
  segs.emplace_back(Point_2(0, 0), Point_2(1-e, 0));
  segs.emplace_back(Point_2(1, 1), Point_2(1, -1));
  test(segs);
  segs.clear();

  // A more complex example
  segs.emplace_back(Point_2(1-e, 1), Point_2(-1-e, -1+2*e));
  segs.emplace_back(Point_2(e/2, e/2), Point_2(1, -1));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(2, 0));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(-2+e, -4));
  segs.emplace_back(Point_2(-2, 2), Point_2(2, 2));
  segs.emplace_back(Point_2(7, 7), Point_2(7+e, 7+e));
  segs.emplace_back(Point_2(5, 7-e), Point_2(9, 7-e));
  test(segs);

  // Minimal intersection subsets of a random examples
  segs.emplace_back(Point_2(0.99999999999992173,-0.9999999999991368),Point_2(0.99999999999987266,-6.016984285966173e-13));
  segs.emplace_back(Point_2(1.000000000000494,-1.0000000000002354),Point_2(0.99999999999950062,8.0391743541535603e-13));
  segs.emplace_back(Point_2(1.0000000000008489,-0.99999999999951794),Point_2(0.99999999999926925,-2.3642280455149055e-13));
  test(segs);
  segs.clear();

  segs.emplace_back(Point_2(1.0000000000008309,9.0061990813884068e-13), Point_2(0.99999999999981259,0.99999999999938394));
  segs.emplace_back(Point_2(1.0000000000005094,-3.1951552372595695e-13), Point_2(1.0000000000004887,1.0000000000006302));
  segs.emplace_back(Point_2(1.0000000000006171,-7.3034227828590228e-13), Point_2(1.0000000000004181,0.99999999999913269));
  segs.emplace_back(Point_2(1.0000000000008491,6.6531536071947998e-14), Point_2(0.9999999999993564,1.0000000000006497));
  segs.emplace_back(Point_2(1.0000000000006859,3.4155608842536406e-13), Point_2(0.99999999999990252,0.99999999999978451));
  segs.emplace_back(Point_2(1.0000000000007121,-1.6363168568197086e-13), Point_2(1.0000000000004778,0.99999999999927969));
  segs.emplace_back(Point_2(1.000000000000753,8.3769047200168284e-13), Point_2(1.0000000000003773,1.0000000000001725));
  segs.emplace_back(Point_2(1.000000000000574,6.6944274898354285e-13), Point_2(1.000000000000427,0.99999999999998845));
  segs.emplace_back(Point_2(1.0000000000006153,-2.798953111315494e-13), Point_2(1.0000000000003044,1.0000000000008289));
  segs.emplace_back(Point_2(1.0000000000008387,5.8214218612051864e-13), Point_2(1.0000000000002918,0.99999999999943212));
  segs.emplace_back(Point_2(1.0000000000007745,-3.9231414307222458e-13), Point_2(1.000000000000294,1.0000000000003408));
  test(segs);
  segs.clear();
}

void test_float_snap_rounding(){
#ifdef BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2
  std::cout << "Test the other API" << std::endl;
#endif

  std::vector< Segment_2 > segs;
  std::vector< Polyline_2 > out;
  FT e(std::pow(2, -60));
  segs.emplace_back(Point_2(1-e, 1), Point_2(-1-e, -1+2*e));
  segs.emplace_back(Point_2(e/2, e/2), Point_2(1, -1));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(2, 0));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(-2+e, -4));
  segs.emplace_back(Point_2(-2, 2), Point_2(2, 2));
  segs.emplace_back(Point_2(7, 7), Point_2(7+e, 7+e));
  segs.emplace_back(Point_2(5, 7-e), Point_2(9, 7-e));

  double_snap_rounding_2(segs.begin(), segs.end(), out);
}

int main(int argc,char *argv[])
{
  CGAL::Random rp;
  CGAL::Random r(argc==1?rp.get_seed():std::stoi(argv[1]));
  std::cout << "random seed = " << r.get_seed() << std::endl;
  std::cout << std::setprecision(17);
  fix_test();
  test_float_snap_rounding();
  test_fully_random(r,1000);
  test_multi_almost_indentical_segments(r,200);
  // test_random_polygons(r,200,10);
  // test_iterative_square_intersection(r,2000);
  return(0);
}
