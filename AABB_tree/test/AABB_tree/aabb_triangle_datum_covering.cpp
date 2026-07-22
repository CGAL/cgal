#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_tree/internal/triangle_datum_covering.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Real_timer.h>

#include <fstream>
#include <iostream>

constexpr int N = 1000;

template <typename Mesh>
void test_no_cover(const Mesh& mesh)
{
  std::cout << "Test WITHOUT cover" << std::endl;

  using GT = typename CGAL::GetGeomTraits<Mesh>::type;
  using FT = typename GT::FT;
  using Point = typename GT::Point_3;
  using Ray_3 = typename GT::Ray_3;
  using Line_3 = typename GT::Line_3;

  using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh, CGAL::Default, CGAL::Tag_true, CGAL::Tag_true>;
  using Traits = CGAL::AABB_traits_3<GT, Primitive>;
  using Tree = CGAL::AABB_tree<Traits>;

  // Build
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);
  std::cout << faces(mesh).size() << " input faces and " << tree.size() << " primitives" << std::endl;

  // Usage
  std::cout << "Traverse..." << std::endl;

  CGAL::Random r; // = CGAL::Random(1337);
  CGAL::Random_points_in_cube_3<Point, CGAL::Creator_uniform_3<FT, Point> > g(1., r);

  std::vector<Point> points;
  points.reserve(N);
  std::copy_n(g, N, std::back_inserter(points));

  CGAL::Real_timer timer;
  timer.start();

  for(int i=0; i<N; ++i)
  {
    bool does_intersect = tree.do_intersect(Line_3(CGAL::ORIGIN, points[i]));
    assert(does_intersect);
  }

  for(int i=0; i<N; ++i)
  {
    auto inter_res = tree.any_intersection(Ray_3(CGAL::ORIGIN, points[i]));
    assert(inter_res);
  }

  for(int i=0; i<N; ++i)
  {
    auto inter_res = tree.first_intersection(Ray_3(CGAL::ORIGIN, points[i]));
    assert(inter_res);
  }

  for(int i=0; i<N; ++i)
  {
    Point res = tree.closest_point(points[i]);
    CGAL_USE(res);
  }

  for(int i=0; i<N; ++i)
    tree.all_intersected_primitives(Line_3(CGAL::ORIGIN, points[i]), CGAL::Emptyset_iterator());

  timer.stop();
  std::cout << "Elapsed: " << timer.time() << " s." << std::endl;
}

template <typename Mesh, typename NamedParameters>
void test_cover(const Mesh& mesh,
                const NamedParameters& np)
{
  std::cout << "Test WITH cover" << std::endl;

  constexpr bool check_with_standard_tree = true;

  using VPM = typename CGAL::GetVertexPointMap<Mesh, NamedParameters>::const_type;
  using GT = typename CGAL::GetGeomTraits<Mesh, NamedParameters>::type;
  using FT = typename GT::FT;
  using Point = typename boost::property_traits<VPM>::value_type;
  using Ray_3 = typename GT::Ray_3;
  using Line_3 = typename GT::Line_3;
  using Triangle_3 = typename GT::Triangle_3;

  using AABB_tree = CGAL::AABB_trees::internal::AABB_covered_triangle_tree<GT, Point>;

  using FG_Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh, CGAL::Default, CGAL::Tag_true, CGAL::Tag_true>;
  using FG_Traits = CGAL::AABB_traits_3<GT, FG_Primitive>;
  using FG_Tree = CGAL::AABB_tree<FG_Traits>;

  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  // Build -----------------------------------------------------------------------------------------
  using CGAL::parameters::get_parameter;
  using CGAL::parameters::choose_parameter;

  VPM vpm = choose_parameter(get_parameter(np, CGAL::internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, mesh));
  GT gt = choose_parameter<GT>(get_parameter(np, CGAL::internal_np::geom_traits));

  AABB_tree tree(diag_length / 100.);
  tree.reserve(num_faces(mesh));

  for(auto f : faces(mesh))
  {
    if(CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(f, mesh, np))
      continue;

    const Point& p0 = get(vpm, target(halfedge(f, mesh), mesh));
    const Point& p1 = get(vpm, target(next(halfedge(f, mesh), mesh), mesh));
    const Point& p2 = get(vpm, source(halfedge(f, mesh), mesh));

    const Triangle_3 tr = gt.construct_triangle_3_object()(p0, p1, p2);
    tree.split_and_insert(tr);
  }

  std::cout << faces(mesh).size() << " input faces and " << tree.size() << " primitives" << std::endl;

  FG_Tree fg_tree(faces(mesh).first, faces(mesh).second, mesh);

  // Usage -----------------------------------------------------------------------------------------
  std::cout << "Traverse..." << std::endl;

  CGAL::Random r; // = CGAL::Random(1337);
  CGAL::Random_points_in_cube_3<Point, CGAL::Creator_uniform_3<FT, Point> > g(1., r);

  std::vector<Point> points;
  points.reserve(N);
  std::copy_n(g, N, std::back_inserter(points));

  using AABB_traits = typename AABB_tree::AABB_traits;
  using AABB_traversal_traits = CGAL::AABB_trees::internal::Covered_tree_traversal_traits<AABB_traits>;

  using Do_intersect_traits = typename AABB_traversal_traits::template Do_intersect_traits<Line_3>;
  using First_intersection_traits = typename AABB_traversal_traits::template First_intersection_traits<Ray_3>;
  using Projection_traits = typename AABB_traversal_traits::Projection_traits;

  using size_type = typename AABB_tree::size_type;
  using Primitive_id = typename AABB_tree::Primitive_id;
  using Counting_iterator = CGAL::internal::AABB_tree::Counting_output_iterator<Primitive_id, size_type>;
  using Listing_primitive_traits = typename AABB_traversal_traits::template Listing_primitive_traits<Line_3, Counting_iterator>;

  CGAL::Real_timer timer;
  timer.start();

  for(int i=0; i<N; ++i)
  {
    Do_intersect_traits do_intersect_traversal_traits(tree.traits());
    Line_3 query(CGAL::ORIGIN, points[i]);
    tree.traversal(query, do_intersect_traversal_traits);
    bool does_intersect = do_intersect_traversal_traits.is_intersection_found();
    assert(does_intersect);

    if(check_with_standard_tree)
      assert(does_intersect == fg_tree.do_intersect(query));
  }

  for(int i=0; i<N; ++i)
  {
    Ray_3 query(CGAL::ORIGIN, points[i]);
    First_intersection_traits first_intersection_traversal_traits(tree.traits());
    tree.traversal(query, first_intersection_traversal_traits);
    auto inter_res = first_intersection_traversal_traits.result();
    assert(inter_res);
  }

  for(int i=0; i<N; ++i)
  {
    Ray_3 query(CGAL::ORIGIN, points[i]);
    auto inter_res = tree.first_intersection(query);
    assert(inter_res);

    if(check_with_standard_tree)
    {
      auto fg_inter_res = fg_tree.first_intersection(query);
      assert(inter_res->first == fg_inter_res->first);
    }
  }

  for(int i=0; i<N; ++i)
  {
    auto hint = tree.best_hint(points[i]);
    Projection_traits projection_traits(hint.first, hint.second, tree.traits());
    Point res = tree.closest_point(points[i]);

    if(check_with_standard_tree)
    {
      Point fg_res = fg_tree.closest_point(points[i]);
      assert(res == fg_res);
    }
  }

  for(int i=0; i<N; ++i)
  {
    Line_3 query(CGAL::ORIGIN, points[i]);

    size_type counter = 0;
    Counting_iterator out(&counter);
    Listing_primitive_traits listing_traversal_traits(out, tree.traits());
    tree.traversal(query, listing_traversal_traits);
    assert(counter > 0);

    if(check_with_standard_tree)
    {
      typename FG_Tree::size_type fg_counter = 0;
      CGAL::internal::AABB_tree::Counting_output_iterator<typename FG_Tree::Primitive_id,
                                                          typename FG_Tree::size_type> fg_out(&fg_counter);
      fg_tree.all_intersected_primitives(query, fg_out);
      assert(counter == fg_counter);
    }
  }

  timer.stop();
  std::cout << "Elapsed: " << timer.time() << " s." << std::endl;
}

template <typename K>
void test(const char* filename)
{
  std::cout << " === Test Kernel " << typeid(K).name() << " === " << std::endl;

  using Point_3 = typename K::Point_3;
  using Mesh = CGAL::Surface_mesh<Point_3>;

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh) ||
     !is_triangle_mesh(mesh))
  {
    std::cerr << "Error: invalid input" << std::endl;
    return;
  }

  test_no_cover(mesh);
  test_cover(mesh, CGAL::parameters::default_values());
}

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  test<CGAL::Exact_predicates_inexact_constructions_kernel>("data/bunny00.off");
  test<CGAL::Exact_predicates_exact_constructions_kernel>("data/bunny00.off");
}

