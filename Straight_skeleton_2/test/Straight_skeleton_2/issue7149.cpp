// #define CGAL_SLS_TEST_ISSUE_7149_DEBUG
#ifdef CGAL_SLS_TEST_ISSUE_7149_DEBUG

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

bool lAppToLog = false ;

void Straight_skeleton_external_trace ( std::string m )
{
  std::ofstream out("sls_log.txt", ( lAppToLog ? std::ios::app | std::ios::ate : std::ios::trunc | std::ios::ate ) );
  out << std::setprecision(17) << m << std::endl << std::flush ;
  lAppToLog = true ;
}
void Straight_skeleton_traits_external_trace ( std::string m )
{
  std::ofstream out("sls_log.txt", ( lAppToLog ? std::ios::app | std::ios::ate : std::ios::trunc | std::ios::ate ) ) ;
  out << std::setprecision(17) << m << std::endl << std::flush ;
  lAppToLog = true ;
}

void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  std::cerr << "CGAL error: " << what << " violation!" << std::endl
       << "Expr: " << expr << std::endl
       << "File: " << file << std::endl
       << "Line: " << line << std::endl;
  if ( msg != nullptr)
      std::cerr << "Explanation:" << msg << std::endl;

  std::exit(1);
}

#define CGAL_STRAIGHT_SKELETON_ENABLE_TRACE 4
#define CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE
#define CGAL_POLYGON_OFFSET_ENABLE_TRACE 4

#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/create_offset_polygons_2.h>

#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/draw_polygon_2.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Random.h>

#include <array>
#include <iostream>

template <typename K, typename PointRange>
void test(const PointRange& points,
          typename K::FT offset)
{
  using FT = typename K::FT;
  using Polygon_2 = CGAL::Polygon_2<K>;

  std::cout << "== Test Kernel: " << typeid(K).name() << ", offset: " << offset << std::endl;

  Polygon_2 pol{std::cbegin(points), std::cend(points)};
  std::cout << "Input polygon is " << (pol.is_simple() ? "simple" : "not simple") << std::endl;

  // For EPICK, construction errors can create polygons with consecutive equal points
  constexpr bool test_output_simplicity =
    (CGAL::is_same_or_derived<CGAL::Field_with_sqrt_tag,
                             typename CGAL::Algebraic_structure_traits<FT>::Algebraic_category>::value &&
     !std::is_floating_point<FT>::value);

  std::vector<Polygon_2> no_holes;
  auto ss_ptr = CGAL::CGAL_SS_i::create_partial_interior_straight_skeleton_2(
                  FT(offset),
                  CGAL::CGAL_SS_i::vertices_begin(pol),
                  CGAL::CGAL_SS_i::vertices_end(pol),
                  no_holes.begin(),
                  no_holes.end(),
                  K());
  assert(ss_ptr);

  std::vector<std::shared_ptr<Polygon_2> > offset_polygons_ptrs =
    CGAL::create_offset_polygons_2<Polygon_2>(FT(offset), CGAL::CGAL_SS_i::dereference(ss_ptr), K());

  std::cout << offset_polygons_ptrs.size() << " polygon(s)" << std::endl;

  if(offset == FT(0.48))
    assert(offset_polygons_ptrs.size() == 1);

  for(const auto& offset_polygon_ptr : offset_polygons_ptrs)
  {
    std::cout << offset_polygon_ptr->size() << " vertices in offset polygon" << std::endl;
    std::cout << "Offset polygon is " << (offset_polygon_ptr->is_simple() ? "simple" : "not simple") << std::endl;
    for(const auto& p : *offset_polygon_ptr)
      std::cout << p << std::endl;

    // CGAL::draw(*offset_polygon_ptr);

    if(test_output_simplicity)
      assert(offset_polygon_ptr->is_simple());
    if(offset == FT(0.48))
      assert(offset_polygon_ptr->size() == 23);
  }
}

template <typename K>
void test(CGAL::Random& r)
{
  using FT = typename K::FT;
  using Polygon_2 = CGAL::Polygon_2<K>;
  using Point_2 = typename K::Point_2;

  // Input
  const std::array<Point_2, 41> points = {{{131.6610, 51.1444},
                                           {132.0460, 50.9782},
                                           {132.0840, 50.9678},
                                           {132.1210, 50.9678},
                                           {132.1480, 50.9574},
                                           {132.2830, 50.9574},
                                           {132.3060, 50.9678},
                                           {132.3800, 50.9678},
                                           {132.6670, 51.0924},
                                           {132.8060, 51.2170},
                                           {132.8040, 51.2274},
                                           {132.8270, 51.2378},
                                           {132.9220, 51.4040},
                                           {132.9940, 51.6324},
                                           {133.0010, 51.6636},
                                           {133.0010, 51.7363},
                                           {133.0080, 51.7674},
                                           {133.0080, 51.8401},
                                           {133.0010, 51.8817},
                                           {133.0010, 51.9544},
                                           {132.9960, 51.9855},
                                           {132.9840, 52.0582},
                                           {132.9770, 52.0998},
                                           {132.9590, 52.1309},
                                           {132.9520, 52.1725},
                                           {132.9300, 52.2348},
                                           {132.9100, 52.2763},
                                           {132.8860, 52.3490},
                                           {132.8680, 52.3802},
                                           {132.7660, 52.5463},
                                           {132.7490, 52.5775},
                                           {132.7210, 52.5982},
                                           {132.7030, 52.6294},
                                           {132.6730, 52.6605},
                                           {132.6280, 52.7125},
                                           {132.5980, 52.7436},
                                           {132.5700, 52.7644},
                                           {132.4830, 52.8371},
                                           {132.4550, 52.8579},
                                           {131.8930, 53.0552},
                                           {131.7120, 53.0344}}};

  Polygon_2 pol(std::cbegin(points), std::cend(points));
  auto ss_ptr = CGAL::create_interior_straight_skeleton_2(pol, K());
  assert(ss_ptr);
  // CGAL::draw(*ss_ptr);

  // get some interesting offset values
  std::set<FT> offsets {{0.48}};

  for(auto it=ss_ptr->vertices_begin(); it!=ss_ptr->vertices_end(); ++it)
  {
    if(it->time() > 0)
    {
      std::cout << "Offset " << it->time() << " at " << it->point() << std::endl;
      offsets.insert(it->time());
    }
  }

  // a couple random values
  for(int i=0; i<10; ++i)
    offsets.insert(FT(r.get_double((std::numeric_limits<double>::min)(),
                      CGAL::to_double(*(offsets.rbegin())))));

  for(FT offset : offsets)
    test<K>(points, offset);
}

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  CGAL::Random r;
  std::cout << "random seed = " << r.get_seed() << std::endl;

#ifdef CGAL_SLS_TEST_ISSUE_7149_DEBUG
  sEnableTrace = true;
#endif

  test<CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt>(r);
  test<CGAL::Exact_predicates_exact_constructions_kernel>(r);
  test<CGAL::Exact_predicates_inexact_constructions_kernel>(r);

  return EXIT_SUCCESS;
}
