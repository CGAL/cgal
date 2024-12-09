#include <iostream>
#include <iomanip>
#include <string>

#define CGAL_SLS_TEST_SPEED_THINGS_UP_FOR_THE_TESTSUITE

//#define CGAL_SLS_PRINT_QUEUE_BEFORE_EACH_POP
//#define CGAL_STRAIGHT_SKELETON_ENABLE_TRACE 100
//#define CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE 10000000
//#define CGAL_STRAIGHT_SKELETON_VALIDITY_ENABLE_TRACE
//#define CGAL_POLYGON_OFFSET_ENABLE_TRACE 10000000

void Straight_skeleton_external_trace(std::string m)
{
  std::cout << std::setprecision(17) << m << std::endl << std::endl ;
}

void Straight_skeleton_traits_external_trace(std::string m)
{
  std::cout << std::setprecision(17) << m << std::endl << std::endl ;
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/draw_polygon_2.h>

#include <memory>

#include <cassert>
#include <iostream>
#include <set>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel            EPECK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  EPECK_w_sqrt;

template <typename K>
void test_API()
{
  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;

  typedef CGAL::Straight_skeleton_2<EPICK>                           Straight_skeleton_EPICK;
  typedef std::shared_ptr<Straight_skeleton_EPICK>                 Straight_skeleton_Ptr_EPICK;

  typedef CGAL::Straight_skeleton_2<K>                               Straight_skeleton;
  typedef std::shared_ptr<Straight_skeleton>                       Straight_skeleton_Ptr;

  Polygon_2 p;
  Straight_skeleton_Ptr_EPICK ss0 = CGAL::create_interior_straight_skeleton_2(p);
  Straight_skeleton_Ptr ss1 = CGAL::create_interior_straight_skeleton_2(p, K());
  Straight_skeleton_Ptr_EPICK ss2 = CGAL::create_exterior_straight_skeleton_2(double(1.01), p);
  Straight_skeleton_Ptr ss3 = CGAL::create_exterior_straight_skeleton_2(int(2), p, K());

  Polygon_with_holes_2 pwh;
  Straight_skeleton_Ptr_EPICK ss4 = CGAL::create_interior_straight_skeleton_2(pwh);
  Straight_skeleton_Ptr ss5 = CGAL::create_interior_straight_skeleton_2(pwh, K());
  Straight_skeleton_Ptr_EPICK ss6 = CGAL::create_exterior_straight_skeleton_2(double(1.01), p);
  Straight_skeleton_Ptr ss7 = CGAL::create_exterior_straight_skeleton_2(int(2), p, K());
}

template <typename K, typename StraightSkeleton>
bool is_valid(const std::shared_ptr<StraightSkeleton>& ss)
{
  typedef typename StraightSkeleton::Traits::Point_2                 Point;
  typedef CGAL::Polygon_2<K>                                         Polygon_2;

  if(!ss)
    return false;

  if(!ss->is_valid())
    return false;

  assert(static_cast<std::size_t>(std::distance(ss->vertices_begin(), ss->vertices_end())) ==
           static_cast<std::size_t>(ss->size_of_vertices()));

  std::set<Point> unique_vertices;
  for(auto vit=ss->vertices_begin(); vit!=ss->vertices_end(); ++vit)
    unique_vertices.insert(vit->point());

  std::cout << unique_vertices.size() << " unique vertices (" << ss->size_of_vertices() << ")" << std::endl;

  // Can't guarantee that the embedding is correct with EPICK or EPECK (non-exact sqrt)
  if(std::is_same<K, EPECK_w_sqrt>::value)
  {
    if(unique_vertices.size() != ss->size_of_vertices())
      return false;

    for(auto hit=ss->halfedges_begin(); hit!=ss->halfedges_end(); ++hit)
      if(hit->vertex()->point() == hit->opposite()->vertex()->point())
        return false;

    for (auto fit=ss->faces_begin(); fit!=ss->faces_end(); ++fit)
    {
      Polygon_2 p;
      auto h = fit->halfedge();
      do
      {
        p.push_back(h->vertex()->point());
        h = h->next();
      }
      while(h != fit->halfedge());

      if(!p.is_simple())
        return false;
    }
  }

  return true;
}

template <typename K>
void test_skeleton(const char* filename,
                   const std::size_t expected_nv = -1,
                   const std::size_t expected_nh = -1,
                   const std::size_t expected_nf = -1)
{
  std::cout << "Construct straight skeleton of input: " << filename << std::endl;
  std::cout << "Kernel: " << typeid(K).name() << std::endl;

  typedef typename K::Point_2                                        Point;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;

  typedef CGAL::Straight_skeleton_2<K>                               Straight_skeleton;
  typedef std::shared_ptr<Straight_skeleton>                       Straight_skeleton_Ptr;

  std::ifstream in(filename);
  assert(in);

  CGAL::IO::set_ascii_mode(in);

  std::vector<Polygon_2> polys;

  int ccb_count = 0;
  in >> ccb_count;
  for(int i=0; i<ccb_count && in; ++i)
  {
    std::vector<Point> poly;

    int v_count = 0;
    in >> v_count;
    for(int j=0; j<v_count && in; ++j)
    {
      double x = 0., y = 0.;
      in >> x >> y;
      poly.emplace_back(x, y);
    }

    if(poly.size() >= 3)
    {
      bool is_simple = CGAL::is_simple_2(poly.begin(), poly.end(), K());
      if(!is_simple)
        std::cerr << "Input polygon not simple (hopefully it is strictly simple...)" << std::endl;

      CGAL::Orientation expected = (i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE);

      const double area = CGAL::to_double(CGAL::polygon_area_2(poly.begin(), poly.end(), K()));
      CGAL::Orientation orientation = area > 0 ? CGAL::COUNTERCLOCKWISE : area < 0 ? CGAL::CLOCKWISE : CGAL::COLLINEAR;

      if(orientation == expected)
        polys.push_back(Polygon_2(poly.begin(), poly.end()));
      else
        polys.push_back(Polygon_2(poly.rbegin(), poly.rend()));
    }
  }

  assert(!polys.empty());

  std::cout <<"Polygon with border of size: " << polys[0].size() << std::endl;
  if(polys.size() > 1)
    std::cout << polys.size() - 1 << " holes" << std::endl;

  Polygon_with_holes_2 p(polys[0]);
  for(std::size_t i=0; i<polys.size()-1; ++i)
    p.add_hole(polys[i+1]);

  std::cout << p.outer_boundary().size() << " vertices" << std::endl;

  Straight_skeleton_Ptr ss = CGAL::create_interior_straight_skeleton_2(p, K());
  assert(ss);
  assert(is_valid<K>(ss));

  std::cout << ss->size_of_vertices() << " vertices" << std::endl;
  std::cout << ss->size_of_halfedges() << " halfedges" << std::endl;
  std::cout << ss->size_of_faces() << " faces" << std::endl;

  if(expected_nv != std::size_t(-1))
  {
    if(expected_nv != ss->size_of_vertices() ||
       expected_nh != ss->size_of_halfedges() ||
       expected_nf != ss->size_of_faces())
    {
      std::cout << "expected: " << expected_nv << " NV " << expected_nh << " NH " << expected_nf << " NF" << std::endl;
      std::cout << "got: " << ss->size_of_vertices() << " NV " << ss->size_of_halfedges() << " NH " << ss->size_of_faces() << " NF" << std::endl;
      assert(false);
    }
  }

  std::cout << "Exterior..." << std::endl;
  Straight_skeleton_Ptr ss_ext = CGAL::create_exterior_straight_skeleton_2(1., p, K());
  assert(ss_ext);
  assert(is_valid<K>(ss_ext));
}

template <typename K>
void test_kernel()
{
//  CGAL_STSKEL_TRAITS_ENABLE_TRACE

#ifndef CGAL_SLS_TEST_SPEED_THINGS_UP_FOR_THE_TESTSUITE
  // test_API<K>();
#endif

  test_skeleton<K>("data/pseudo_split_0.poly", 13, 40, 8);
  test_skeleton<K>("data/pseudo_split_1.poly", 21, 68, 12);
  test_skeleton<K>("data/pseudo_split_2.poly", 21, 68, 12);
  test_skeleton<K>("data/pseudo_split_3.poly", 25, 86, 16);
  test_skeleton<K>("data/pseudo_split_4.poly", 25, 86, 16);
  test_skeleton<K>("data/pseudo_split_5.poly"/*, 29, 104, 20*/); // almost duplicates
  test_skeleton<K>("data/pseudo_split_5b.poly"/*, 29, 104, 20*/); // almost duplicates
  test_skeleton<K>("data/pseudo_split_6.poly", 13, 40, 8);
  test_skeleton<K>("data/pseudo_split_7.poly", 65, 216, 44);
  test_skeleton<K>("data/pseudo_split_8.poly");
  test_skeleton<K>("data/pseudo_split_9.poly");
  test_skeleton<K>("data/pseudo_split_10.poly", 13, 40, 8);
  test_skeleton<K>("data/pseudo_split_11.poly", 14, 44, 9);
  test_skeleton<K>("data/pseudo_split_12.poly");
  test_skeleton<K>("data/pseudo_split_13b.poly", 11, 34, 7);
  test_skeleton<K>("data/pseudo_split_13.poly", 11, 34, 7);

  test_skeleton<K>("data/1_Example.poly", 119, 364, 64);
  test_skeleton<K>("data/1_Example_Working.poly", 104, 318, 56);
  test_skeleton<K>("data/2_Example.poly", 37, 112, 20);
  test_skeleton<K>("data/5-SPOKE2.poly"/*, 32, 102, 22*/); // weird almost-duplicates
  test_skeleton<K>("data/5-SPOKE.poly"/*, 32, 102, 22*/);
  test_skeleton<K>("data/7-SPOKE.poly");
  test_skeleton<K>("data/alley_0.poly", 18, 58, 12);
  test_skeleton<K>("data/alley_1.poly", 20, 62, 12);
  test_skeleton<K>("data/alley_2.poly", 26, 78, 12);
  test_skeleton<K>("data/alley_3.poly", 34, 102, 16);
  test_skeleton<K>("data/AlmostClosed.poly");
  test_skeleton<K>("data/A.poly", 58, 174, 29);
  test_skeleton<K>("data/closer_edge_event_0.poly", 11, 34, 7);
  test_skeleton<K>("data/closer_edge_event_1.poly", 11, 34, 7);
  test_skeleton<K>("data/consecutive_coincident_vertices_0.poly", 6, 18, 4);
  test_skeleton<K>("data/consecutive_coincident_vertices_1.poly", 6, 18, 4);
  test_skeleton<K>("data/consecutive_coincident_vertices_2.poly", 6, 18, 4);
  test_skeleton<K>("data/consecutive_coincident_vertices_3.poly", 6, 18, 4);
  test_skeleton<K>("data/consecutive_coincident_vertices_4.poly", 6, 18, 4);
  test_skeleton<K>("data/degenerate0.poly", 8, 24, 5);
  test_skeleton<K>("data/degenerate0a.poly", 8, 24, 5);
  test_skeleton<K>("data/degenerate1.poly", 9, 28, 6);
  test_skeleton<K>("data/degenerate2.poly", 7, 22, 5);
  test_skeleton<K>("data/degenerate3.poly", 10, 32, 7);
  test_skeleton<K>("data/degenerate4.poly", 11, 36, 8);
  test_skeleton<K>("data/degenerate5a.poly", 10, 30, 6);
  test_skeleton<K>("data/degenerate5.poly", 10, 30, 6);
  test_skeleton<K>("data/degenerate6.poly", 8, 26, 6);
  test_skeleton<K>("data/degenerate7.poly", 12, 36, 7);
  test_skeleton<K>("data/degenerate8.poly", 15, 52, 12);
  test_skeleton<K>("data/degenerate9.poly", 16, 48, 9);
  test_skeleton<K>("data/degenerate10.poly", 19, 58, 11);
  test_skeleton<K>("data/degenerate11.poly", 17, 56, 12);
  test_skeleton<K>("data/degenerate12.poly", 17, 52, 10);
  test_skeleton<K>("data/degenerate13.poly", 9, 32, 8);
  test_skeleton<K>("data/degenerate20.poly");
  test_skeleton<K>("data/degenerate21.poly");
  test_skeleton<K>("data/degenerate22.poly", 36, 108, 18);
  test_skeleton<K>("data/degenerate22b.poly", 35, 106, 18);
  test_skeleton<K>("data/degenerate22c.poly", 29, 88, 15);
  test_skeleton<K>("data/degenerate24.poly", 36, 108, 18);
  test_skeleton<K>("data/degenerate25.poly", 19, 58, 10);
  test_skeleton<K>("data/degenerate26.poly", 24, 74, 14);
  test_skeleton<K>("data/degenerate27b.poly", 19, 58, 11);
  test_skeleton<K>("data/degenerate27c.poly", 19, 58, 11);
  test_skeleton<K>("data/degenerate27d.poly", 25, 76, 14);
  test_skeleton<K>("data/degenerate27e.poly", 26, 80, 15);
  test_skeleton<K>("data/degenerate27.poly", 18, 56, 11);
  test_skeleton<K>("data/degenerate28aa.poly"/*, 58, 178, 32*/); // should be 58, but almost-duplicate so 59
  test_skeleton<K>("data/degenerate28a.poly"/*, 58, 178, 32*/); // same as above
  test_skeleton<K>("data/degenerate28bb.poly");
  test_skeleton<K>("data/degenerate28b.poly");
  test_skeleton<K>("data/degenerate28c.poly");
  test_skeleton<K>("data/degenerate28x.poly");
  test_skeleton<K>("data/degenerate29.poly");
  test_skeleton<K>("data/degenerate_multinode0.poly", 41, 136, 28);
  test_skeleton<K>("data/double_edge.poly", 7, 22, 5);
  test_skeleton<K>("data/double_edge_0.poly", 7, 22, 5);
  test_skeleton<K>("data/double_edge_1.poly", 7, 22, 5);
  test_skeleton<K>("data/double_edge_2.poly", 35, 108, 19);
  test_skeleton<K>("data/double_split.poly", 12, 36, 7);
  test_skeleton<K>("data/equal_times_0.poly", 16, 48, 9);
  test_skeleton<K>("data/ExtraEdge_1.poly");
  test_skeleton<K>("data/ExtraEdge_2.poly"/*, 41, 130, 25*/); // almost duplicates
  test_skeleton<K>("data/hole.poly", 16, 48, 8);
  test_skeleton<K>("data/inputcircle.poly"/* 61, 188, 34*/); // almost duplicates
  test_skeleton<K>("data/inputsquare.poly", 29, 104, 20);
  test_skeleton<K>("data/inputsquare2.poly", 21, 68, 12);
  test_skeleton<K>("data/many_holes.poly");
  test_skeleton<K>("data/masked_double_split.poly", 10, 30, 6);
  test_skeleton<K>("data/multinode0.poly", 29, 96, 20);
  test_skeleton<K>("data/multinode1.poly", 29, 96, 20);
  test_skeleton<K>("data/near_degenerate_0.poly", 10, 30, 6);
  test_skeleton<K>("data/near_degenerate_1.poly", 14, 42, 8);
  test_skeleton<K>("data/nearly_collinear.poly");
  test_skeleton<K>("data/parallels0_b.poly");
  test_skeleton<K>("data/parallels0.poly");
  test_skeleton<K>("data/parallels_1.poly", 30, 106, 20);
  test_skeleton<K>("data/poly6.poly"/*, 15, 48, 10*/); // almost duplicates
  test_skeleton<K>("data/rect_4_spokes.poly", 21, 72, 16);
  test_skeleton<K>("data/rectangle.poly", 6, 18, 4);
  test_skeleton<K>("data/region_4.poly", 4, 12, 3);
  test_skeleton<K>("data/rombus_4_spokes.poly", 17, 56, 12);
  test_skeleton<K>("data/issue5177.poly");

  // The embedding of those below is bad when using EPICK
  test_skeleton<K>("data/poly4.poly"/*, 15, 48, 10*/); // almost duplicates
  test_skeleton<K>("data/poly4b.poly"/*, 15, 48, 10*/); // almost duplicates
  test_skeleton<K>("data/Detmier.poly");
  test_skeleton<K>("data/Detmier_b.poly"/*, 28, 86, 16*/); // almost duplicates
  test_skeleton<K>("data/Detmier_c.poly");
  test_skeleton<K>("data/Detmier_d.poly"/*, 27, 84, 16*/); // almost duplicates
  test_skeleton<K>("data/Detmier_e.poly");

#ifndef CGAL_SLS_TEST_SPEED_THINGS_UP_FOR_THE_TESTSUITE
  test_skeleton<K>("data/sample_0.poly");
  test_skeleton<K>("data/sample_101.poly");
  test_skeleton<K>("data/sample_102.poly");
  test_skeleton<K>("data/sample_147.poly");
  test_skeleton<K>("data/sample_1.poly", 24, 72, 12);
  test_skeleton<K>("data/sample_235.poly");
  test_skeleton<K>("data/sample_298.poly");
  test_skeleton<K>("data/sample_2.poly");
  test_skeleton<K>("data/sample2.poly", 12, 36, 7);
  test_skeleton<K>("data/sample_319.poly");
  test_skeleton<K>("data/sample_325.poly");
  test_skeleton<K>("data/sample_333.poly");
  test_skeleton<K>("data/sample_3.poly", 16, 48, 8);
  test_skeleton<K>("data/sample3.poly");
  test_skeleton<K>("data/sample_46.poly");
  test_skeleton<K>("data/sample_4.poly");
  test_skeleton<K>("data/sample_5.poly");
  test_skeleton<K>("data/sample_638.poly");
  test_skeleton<K>("data/sample_698.poly");
  test_skeleton<K>("data/sample_6.poly", 18, 54, 10);
  test_skeleton<K>("data/sample_73.poly");
  test_skeleton<K>("data/sample_85.poly");
  test_skeleton<K>("data/sample.poly", 22, 66, 12);
  test_skeleton<K>("data/simple_0.poly");
  test_skeleton<K>("data/simple_1.poly", 24, 72, 12);
  test_skeleton<K>("data/simple_2.poly");
  test_skeleton<K>("data/simple_3.poly");
  test_skeleton<K>("data/single_split.poly", 8, 24, 5);
  test_skeleton<K>("data/split_at_end_0.poly", 11, 34, 7);
  test_skeleton<K>("data/split_at_end_1.poly", 11, 36, 8);
  test_skeleton<K>("data/split_at_end_2.poly", 11, 34, 7);
  test_skeleton<K>("data/split_at_zero_0.poly", 20, 60, 11);
  test_skeleton<K>("data/square.poly", 5, 16, 4);
  test_skeleton<K>("data/triangle.poly", 4, 12, 3);
  test_skeleton<K>("data/star.poly", 9, 32, 8);
  test_skeleton<K>("data/StrayCenterlines.poly");
  test_skeleton<K>("data/wheel_13_spokes.poly");
  test_skeleton<K>("data/wheel_14_spokes.poly");
  test_skeleton<K>("data/wheel_15_spokes.poly");
  test_skeleton<K>("data/wheel_16_spokes.poly");
  test_skeleton<K>("data/wheel_16_spokes_b.poly");
  test_skeleton<K>("data/wheel_128_spokes.poly");
  test_skeleton<K>("data/wiggly_03_cgal.poly");
  test_skeleton<K>("data/WingChiu.poly");

  test_skeleton<K>("data/large_1.poly");
  test_skeleton<K>("data/large_2.poly");
  test_skeleton<K>("data/large_3.poly");
  test_skeleton<K>("data/large_4.poly");
  test_skeleton<K>("data/complex_0.poly");
  test_skeleton<K>("data/complex_1.poly");
  test_skeleton<K>("data/complex_2.poly");
  test_skeleton<K>("data/complex_3.poly");
  test_skeleton<K>("data/complex_4.poly");
  test_skeleton<K>("data/complex_5.poly");
  test_skeleton<K>("data/inputc.poly");
  test_skeleton<K>("data/inputd1.poly");
  test_skeleton<K>("data/inputd.poly");
  test_skeleton<K>("data/inputG.poly");
  test_skeleton<K>("data/input_K.poly");
  test_skeleton<K>("data/inputPa.poly");
  test_skeleton<K>("data/inputP.poly");
  test_skeleton<K>("data/inputq.poly");
  test_skeleton<K>("data/inputT.poly");
  test_skeleton<K>("data/inputu.poly");
  test_skeleton<K>("data/inputq1.poly", 21, 64, 12);
#endif
}

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  test_kernel<EPICK>();
  test_kernel<EPECK>();
  test_kernel<EPECK_w_sqrt>();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
