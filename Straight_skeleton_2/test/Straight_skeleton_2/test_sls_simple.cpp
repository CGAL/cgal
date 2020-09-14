#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include "print.h"

#include <CGAL/Polygon_2.h>
#include <CGAL/draw_polygon_2.h>

#include <boost/shared_ptr.hpp>

#include <cassert>
#include <iostream>
#include <set>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel            EPECK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  EPECK_w_sqrt;

template <typename StraightSkeleton>
bool is_valid(const boost::shared_ptr<StraightSkeleton>& ss)
{
  typedef typename StraightSkeleton::Traits::Point_2 Point;

  if(!ss->is_valid())
    return false;

  std::set<Point> unique_vertices;
  for(auto vit=ss->vertices_begin(); vit!=ss->vertices_end(); ++vit)
    unique_vertices.insert(vit->point());
  if(unique_vertices.size() != ss->size_of_vertices())
    return false;

  for(auto hit=ss->halfedges_begin(); hit!=ss->halfedges_end(); ++hit)
    if(hit->vertex()->point() == hit->opposite()->vertex()->point())
      return false;

  return true;
}

template <typename K>
void test_skeleton(const char* filename,
                   const std::size_t expected_nv = -1,
                   const std::size_t expected_nh = -1,
                   const std::size_t expected_nf = -1)
{
  std::cout << "Construct straight skeleton of input: " << filename << std::endl;

  typedef typename K::FT                                             FT;
  typedef typename K::Point_2                                        Point;

  typedef CGAL::Polygon_2<K>                                         Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>                              Polygon_with_holes_2;
  typedef boost::shared_ptr<Polygon_with_holes_2>                    Polygon_with_holes_2_ptr;
  typedef std::vector<Polygon_with_holes_2_ptr>                      Polygon_with_holes_2_ptr_container;

  typedef CGAL::Straight_skeleton_2<K>                               Straight_skeleton;
  typedef boost::shared_ptr<Straight_skeleton>                       Straight_skeleton_Ptr;

  std::ifstream in(filename);
  assert(in);

  CGAL::set_ascii_mode(in);

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
      poly.push_back(Point(x, y));
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

  Polygon_with_holes_2 p(polys[0]);
  for(std::size_t i=0; i<polys.size()-1; ++i)
    p.add_hole(polys[1]);

  Straight_skeleton_Ptr ss = CGAL::create_interior_straight_skeleton_2(p);

  std::cout << ss->size_of_vertices() << " vertices" << std::endl;
  std::cout << ss->size_of_halfedges() << " halfedges" << std::endl;
  std::cout << ss->size_of_faces() << " faces" << std::endl;

  if(expected_nv != std::size_t(-1))
    assert(expected_nv == ss->size_of_vertices());
  if(expected_nh != std::size_t(-1))
    assert(expected_nh == ss->size_of_halfedges());
  if(expected_nf != std::size_t(-1))
    assert(expected_nf == ss->size_of_faces());

  assert(is_valid(ss));

  ss = CGAL::create_interior_straight_skeleton_2(p, K());

  if(expected_nv != std::size_t(-1))
    assert(expected_nv == ss->size_of_vertices());
  if(expected_nh != std::size_t(-1))
    assert(expected_nh == ss->size_of_halfedges());
  if(expected_nf != std::size_t(-1))
    assert(expected_nf == ss->size_of_faces());

  assert(is_valid(ss));


  // @todo test API range
}

template <typename K>
void test_kernel()
{
  std::cout.precision(17);
  std::cerr.precision(17);

  // Disk, NV + NF - 1 = NE
  test_skeleton<K>("data/1_Example.poly", 119, 364, 64);
  test_skeleton<K>("data/1_Example_Working.poly", 104, 318, 56);
  test_skeleton<K>("data/2_Example.poly", 37, 112, 20);
  test_skeleton<K>("data/5-SPOKE2.poly"/*, 32, 102, 22*/); // weird almost-duplicates
  test_skeleton<K>("data/5-SPOKE.poly"/*, 32, 102, 22*/);
  test_skeleton<K>("data/7-SPOKE.poly"); // @todo recompute
  test_skeleton<K>("data/alley_0.poly", 18, 58, 12);
  test_skeleton<K>("data/alley_1.poly", 20, 62, 12);
  test_skeleton<K>("data/alley_2.poly");
  test_skeleton<K>("data/alley_3.poly");
  test_skeleton<K>("data/AlmostClosed.poly");
  test_skeleton<K>("data/A.poly");
  test_skeleton<K>("data/closer_edge_event_0.poly");
  test_skeleton<K>("data/closer_edge_event_1.poly");
  test_skeleton<K>("data/complex_0.poly");
  test_skeleton<K>("data/complex_1.poly");
  test_skeleton<K>("data/complex_2.poly");
  test_skeleton<K>("data/complex_3.poly");
  test_skeleton<K>("data/complex_4.poly");
  test_skeleton<K>("data/complex_5.poly");
  test_skeleton<K>("data/consecutive_coincident_vertices_0.poly");
  test_skeleton<K>("data/consecutive_coincident_vertices_1.poly");
  test_skeleton<K>("data/consecutive_coincident_vertices_2.poly");
  test_skeleton<K>("data/consecutive_coincident_vertices_3.poly");
  test_skeleton<K>("data/consecutive_coincident_vertices_4.poly");
  test_skeleton<K>("data/degenerate0a.poly");
  test_skeleton<K>("data/degenerate0.poly");
  test_skeleton<K>("data/degenerate10.poly");
  test_skeleton<K>("data/degenerate11.poly");
  test_skeleton<K>("data/degenerate12.poly");
  test_skeleton<K>("data/degenerate13.poly");
  test_skeleton<K>("data/degenerate1.poly");
  test_skeleton<K>("data/degenerate20.poly");
  test_skeleton<K>("data/degenerate21.poly");
  test_skeleton<K>("data/degenerate22b.poly");
  test_skeleton<K>("data/degenerate22c.poly");
  test_skeleton<K>("data/degenerate22.poly");
  test_skeleton<K>("data/degenerate24.poly");
  test_skeleton<K>("data/degenerate25.poly");
  test_skeleton<K>("data/degenerate26.poly");
  test_skeleton<K>("data/degenerate27b.poly");
  test_skeleton<K>("data/degenerate27c.poly");
  test_skeleton<K>("data/degenerate27d.poly");
  test_skeleton<K>("data/degenerate27e.poly");
  test_skeleton<K>("data/degenerate27.poly");
  test_skeleton<K>("data/degenerate28aa.poly");
  test_skeleton<K>("data/degenerate28a.poly");
  test_skeleton<K>("data/degenerate28b.poly");
  test_skeleton<K>("data/degenerate28c.poly");
  test_skeleton<K>("data/degenerate28x.poly");
  test_skeleton<K>("data/degenerate2.poly");
  test_skeleton<K>("data/degenerate3.poly");
  test_skeleton<K>("data/degenerate4.poly");
  test_skeleton<K>("data/degenerate5a.poly");
  test_skeleton<K>("data/degenerate5.poly");
  test_skeleton<K>("data/degenerate6.poly");
  test_skeleton<K>("data/degenerate7.poly");
  test_skeleton<K>("data/degenerate8.poly");
  test_skeleton<K>("data/degenerate9.poly");
  test_skeleton<K>("data/degenerate_multinode0.poly");
  test_skeleton<K>("data/Detmier_b.poly");
  test_skeleton<K>("data/Detmier_c.poly");
  test_skeleton<K>("data/Detmier_d.poly");
  test_skeleton<K>("data/Detmier_e.poly");
  test_skeleton<K>("data/Detmier.poly");
  test_skeleton<K>("data/double_edge_0.poly");
  test_skeleton<K>("data/double_edge_1.poly");
  test_skeleton<K>("data/double_edge_2.poly");
  test_skeleton<K>("data/double_edge.poly");
  test_skeleton<K>("data/double_split.poly");
  test_skeleton<K>("data/equal_times_0.poly");
  test_skeleton<K>("data/ExtraEdge_1.poly");
  test_skeleton<K>("data/ExtraEdge_2.poly");
  test_skeleton<K>("data/hole.poly");
  test_skeleton<K>("data/inputcircle.poly");
  test_skeleton<K>("data/inputc.poly");
  test_skeleton<K>("data/inputd1.poly");
  test_skeleton<K>("data/inputd.poly");
  test_skeleton<K>("data/inputG.poly");
  test_skeleton<K>("data/input_K.poly");
  test_skeleton<K>("data/inputPa.poly");
  test_skeleton<K>("data/inputP.poly");
  test_skeleton<K>("data/inputq1.poly");
  test_skeleton<K>("data/inputq.poly");
  test_skeleton<K>("data/inputsquare2.poly");
  test_skeleton<K>("data/inputsquare.poly");
  test_skeleton<K>("data/inputT.poly");
  test_skeleton<K>("data/inputu.poly");
  test_skeleton<K>("data/issue3382_bis.txt");
  test_skeleton<K>("data/issue3382_ter.txt");
  test_skeleton<K>("data/large_1.poly");
  test_skeleton<K>("data/large_2.poly");
  test_skeleton<K>("data/large_3.poly");
  test_skeleton<K>("data/large_4.poly");
  test_skeleton<K>("data/many_holes.poly");
  test_skeleton<K>("data/masked_double_split.poly");
  test_skeleton<K>("data/multinode0.poly");
  test_skeleton<K>("data/multinode1.poly");
  test_skeleton<K>("data/near_degenerate_0.poly");
  test_skeleton<K>("data/near_degenerate_1.poly");
  test_skeleton<K>("data/nearly_collinear.poly");
  test_skeleton<K>("data/parallels0_b.poly");
  test_skeleton<K>("data/parallels0.poly");
  test_skeleton<K>("data/parallels_1.poly");
  test_skeleton<K>("data/poly4b.poly");
  test_skeleton<K>("data/poly4.poly");
  test_skeleton<K>("data/poly6.poly");
  test_skeleton<K>("data/pseudo_split_0.poly");
  test_skeleton<K>("data/pseudo_split_10.poly");
  test_skeleton<K>("data/pseudo_split_11.poly");
  test_skeleton<K>("data/pseudo_split_12.poly");
  test_skeleton<K>("data/pseudo_split_13b.poly");
  test_skeleton<K>("data/pseudo_split_13.poly");
  test_skeleton<K>("data/pseudo_split_1.poly");
  test_skeleton<K>("data/pseudo_split_2.poly");
  test_skeleton<K>("data/pseudo_split_3.poly");
  test_skeleton<K>("data/pseudo_split_4.poly");
  test_skeleton<K>("data/pseudo_split_5b.poly");
  test_skeleton<K>("data/pseudo_split_5.poly");
  test_skeleton<K>("data/pseudo_split_6.poly");
  test_skeleton<K>("data/pseudo_split_7.poly");
  test_skeleton<K>("data/pseudo_split_8.poly");
  test_skeleton<K>("data/pseudo_split_9.poly");
  test_skeleton<K>("data/rect_4_spokes.poly");
  test_skeleton<K>("data/rectangle.poly");
  test_skeleton<K>("data/region_4.poly");
  test_skeleton<K>("data/rombus_4_spokes.poly");
  test_skeleton<K>("data/sample_0.poly");
  test_skeleton<K>("data/sample_101.poly");
  test_skeleton<K>("data/sample_102.poly");
  test_skeleton<K>("data/sample_147.poly");
  test_skeleton<K>("data/sample_1.poly");
  test_skeleton<K>("data/sample_235.poly");
  test_skeleton<K>("data/sample_298.poly");
  test_skeleton<K>("data/sample_2.poly");
  test_skeleton<K>("data/sample2.poly");
  test_skeleton<K>("data/sample_319.poly");
  test_skeleton<K>("data/sample_325.poly");
  test_skeleton<K>("data/sample_333.poly");
  test_skeleton<K>("data/sample_3.poly");
  test_skeleton<K>("data/sample3.poly");
  test_skeleton<K>("data/sample_46.poly");
  test_skeleton<K>("data/sample_4.poly");
  test_skeleton<K>("data/sample_5.poly");
  test_skeleton<K>("data/sample_638.poly");
  test_skeleton<K>("data/sample_698.poly");
  test_skeleton<K>("data/sample_6.poly");
  test_skeleton<K>("data/sample_73.poly");
  test_skeleton<K>("data/sample_85.poly");
  test_skeleton<K>("data/sample.poly");
  test_skeleton<K>("data/simple_0.poly");
  test_skeleton<K>("data/simple_1.poly");
  test_skeleton<K>("data/simple_2.poly");
  test_skeleton<K>("data/simple_3.poly");
  test_skeleton<K>("data/single_split.poly");
  test_skeleton<K>("data/split_at_end_0.poly");
  test_skeleton<K>("data/split_at_end_1.poly");
  test_skeleton<K>("data/split_at_end_2.poly");
  test_skeleton<K>("data/split_at_zero_0.poly");
  test_skeleton<K>("data/square.poly");
  test_skeleton<K>("data/star.poly");
  test_skeleton<K>("data/StrayCenterlines.poly");
  test_skeleton<K>("data/triangle.poly");
  test_skeleton<K>("data/wheel_128_spokes.poly");
  test_skeleton<K>("data/wheel_13_spokes.poly");
  test_skeleton<K>("data/wheel_14_spokes.poly");
  test_skeleton<K>("data/wheel_15_spokes.poly");
  test_skeleton<K>("data/wheel_16_spokes_b.poly");
  test_skeleton<K>("data/wheel_16_spokes.poly");
  test_skeleton<K>("data/wiggly_03_cgal.poly");
  test_skeleton<K>("data/WingChiu.poly");
}

int main(int, char**)
{
  test_kernel<EPICK>();
  test_kernel<EPECK>();
  test_kernel<EPECK_w_sqrt>();
}
