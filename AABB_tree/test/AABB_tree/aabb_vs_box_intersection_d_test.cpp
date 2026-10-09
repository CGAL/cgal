#include <CGAL/box_intersection_d.h>
#include <CGAL/AABB_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_trees/intersection.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <iostream>
#include <set>
#include <cstdio>
#include <cassert>

#include <CGAL/Random.h>

using Box = CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::size_t>;
using BoxRange = std::vector<Box>;
using K = CGAL::Simple_cartesian<double>;

template <class Iterator>
struct Iterator_to_bbox_property_map{
  typedef Iterator key_type;
  typedef CGAL::Bbox_3 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  typedef Iterator_to_bbox_property_map<Iterator> Self;

  inline friend reference
  get(Self, key_type it)
  {
    return it->bbox();
  }
};

template <class GeomTraits, class Iterator>
struct Point_in_bbox_3_iterator_property_map{
  typedef Iterator key_type;
  typedef typename GeomTraits::Point_3 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  typedef Point_in_bbox_3_iterator_property_map<GeomTraits, Iterator> Self;

  inline friend reference
  get(Self, key_type it)
  {
    return typename GeomTraits::Point_3( it->bbox().xmin(), it->bbox().ymin(), it->bbox().zmin() );
  }
};

template < class GeomTraits, class Iterator>
struct AABB_box_primitive_3 : public CGAL::AABB_primitive< Iterator,
                                                           Iterator_to_bbox_property_map<Iterator>,
                                                           Point_in_bbox_3_iterator_property_map<GeomTraits, Iterator>,
                                                           CGAL::Tag_false, CGAL::Tag_false >
{
  typedef CGAL::AABB_primitive< Iterator,
                                Iterator_to_bbox_property_map<Iterator>,
                                Point_in_bbox_3_iterator_property_map<GeomTraits, Iterator>,
                                CGAL::Tag_false, CGAL::Tag_false > Base;
  AABB_box_primitive_3(Iterator it) : Base(it){}
};
using Primitive = AABB_box_primitive_3<K, BoxRange::const_iterator>;
using AABB_traits = CGAL::AABB_traits_3<K, Primitive>;
using Tree = CGAL::AABB_tree<AABB_traits>;

void test(BoxRange& r1, BoxRange& r2)
{
  int aabb_count = 0,
      box_intersection_count = 0;

  Tree tree1(r1.begin(), r1.end());
  Tree tree2(r2.begin(), r2.end());

  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, boost::make_function_output_iterator( [&](auto&&) { ++aabb_count; }) );
  CGAL::box_intersection_d(r1.begin(), r1.end(), r2.begin(), r2.end(), [&](auto&, auto&){ ++box_intersection_count; });

  assert(aabb_count == box_intersection_count);
}

void random_test(std::size_t n = 1000)
{
  CGAL::Random r;

  BoxRange boxes1, boxes2;
  boxes1.reserve(n);
  boxes2.reserve(n);

  auto generate_boxes = [&](BoxRange& boxes)
  {
    for(std::size_t i = 0; i < n; ++i){
      const double x1 = r.get_double(0,10),
                   y1 = r.get_double(0,10),
                   z1 = r.get_double(0,10);
      const double x2 = r.get_double(0,10),
                   y2 = r.get_double(0,10),
                   z2 = r.get_double(0,10);
      boxes.emplace_back(CGAL::Bbox_3((std::min)(x1,x2), (std::min)(y1,y2), (std::min)(z1,z2),
                                      (std::max)(x1,x2), (std::max)(y1,y2), (std::max)(z1,z2)), i);
    }
  };
  auto generate_unit_boxes = [&](BoxRange& boxes)
  {
    for(std::size_t i = 0; i < n; ++i){
      const double x = r.get_double(0,9),
                   y = r.get_double(0,9),
                   z = r.get_double(0,9);
      boxes.emplace_back(CGAL::Bbox_3(x, y, z, x+1, y+1, z+1), i);
    }
  };
  generate_boxes(boxes1);
  generate_boxes(boxes2);
  test(boxes1, boxes2);
  boxes1.clear(); boxes2.clear();
  generate_unit_boxes(boxes1);
  generate_unit_boxes(boxes2);
  test(boxes1, boxes2);
}

template <typename Path>
void test_from_data(Path filename1, Path filename2)
{
  using Mesh = CGAL::Surface_mesh<K::Point_3>;

  Mesh mesh1, mesh2;
  if(!CGAL::IO::read_polygon_mesh(filename1, mesh1) || !CGAL::IO::read_polygon_mesh(filename2, mesh2)){
    std::cout << "Files not found." << std::endl;
    return;
  }

  auto make_boxes = [](const Mesh& mesh){
    BoxRange boxes;
    std::size_t i = 0;
    for(auto f: mesh.faces())
      boxes.emplace_back(CGAL::Polygon_mesh_processing::face_bbox(f, mesh), i++);
    return boxes;
  };

  BoxRange boxes1 = make_boxes(mesh1);
  BoxRange boxes2 = make_boxes(mesh2);
  test(boxes1, boxes2);
}

int main()
{
  test_from_data(CGAL::data_file_path("meshes/diplodocus.off"), CGAL::data_file_path("meshes/bunny00.off"));
  test_from_data(CGAL::data_file_path("meshes/bear.off"), CGAL::data_file_path("meshes/man.off"));
  test_from_data(CGAL::data_file_path("meshes/sphere.off"), CGAL::data_file_path("meshes/elephant.off"));
  test_from_data(CGAL::data_file_path("meshes/cow.off"), CGAL::data_file_path("meshes/cross.off"));
  random_test();
  return 0;
}