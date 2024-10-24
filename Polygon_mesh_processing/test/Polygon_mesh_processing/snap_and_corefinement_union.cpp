#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/internal/surface_snapping.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Polygon_mesh.h>
#include <CGAL/Shape_regularization/regularize_planes.h>

#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Kernel_traits.h>

#include <CGAL/IO/write_ply_points.h>
#include <CGAL/IO/Color.h>

#include <CGAL/Index_property_map.h>

#include <CGAL/Kernel_traits.h>

#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

template <class Point_3, class Iterator, class Mesh, class VPM>
struct Point_from_halfedge_map {
  //classical typedefs
  using key_type = Iterator;
  using value_type = Point_3;
  using reference = Point_3&;
  using category = boost::readable_property_map_tag;
  using Self = Point_from_halfedge_map<Point_3, Iterator, Mesh, VPM>;

  Point_from_halfedge_map() : vpms(nullptr), meshes(nullptr) { CGAL_assertion(false); }

  Point_from_halfedge_map(const std::vector<VPM>* vpms, const std::vector<Mesh*>* meshes) : vpms(vpms), meshes(meshes) {}

  inline friend reference
    get(Self s, key_type it) {
    return get((*s.vpms)[it.first], target(it.second, *(*s.meshes)[it.first]));
  }

  const std::vector<VPM>* vpms;
  const std::vector<Mesh*>* meshes;
  //std::vector<std::reference_wrapper<Mesh>>
};

template <class Primitive>
struct Plane_map {
  //classical typedefs
  typedef Primitive key_type;
  typedef typename Primitive::first_type value_type;
  typedef typename Primitive::first_type& reference;
  typedef boost::readable_property_map_tag category;
  typedef Plane_map<Primitive> Self;

  inline friend reference get(Self, key_type it)
  {
    return it.first;
  }

  inline friend reference put(Self, key_type it, value_type value)
  {
    return it.first = value;
  }
};

template<typename PrimitiveRange>
std::vector<K::Plane_3> unique_planes(const PrimitiveRange& prims) {
  std::vector<K::Plane_3> planes;
  for (const auto& p : prims) {
    bool exists = false;
    for (std::size_t i = 0; i < planes.size(); i++) {
      if (planes[i] == p || planes[i].opposite() == p) {
        exists = true;
        break;
      }
    }
    if (!exists)
      planes.push_back(p);
  }

  return planes;
}

template<typename TriangleMesh1, typename TriangleMesh2, typename FT,
  typename NamedParameters1 = CGAL::parameters::Default_named_parameters,
  typename NamedParameters2 = CGAL::parameters::Default_named_parameters>
void coregularize(TriangleMesh1 &mesh1, TriangleMesh2 &mesh2, FT epsilon, NamedParameters1 np1 = CGAL::parameters::default_values(), NamedParameters2 np2 = CGAL::parameters::default_values()) {
  static_assert(std::is_same<typename CGAL::property_map_selector<TriangleMesh1, boost::vertex_point_t>::type, typename CGAL::property_map_selector<TriangleMesh2, boost::vertex_point_t>::type>::value);
  using VPM = typename CGAL::property_map_selector<TriangleMesh1, boost::vertex_point_t>::type;
  //using VPM2 = typename CGAL::property_map_selector<TriangleMesh2, boost::vertex_point_t>::type;

  static_assert(std::is_same<TriangleMesh1, TriangleMesh2>::value);
  using Mesh = TriangleMesh1;

  using K = typename CGAL::Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;
  using Point_3 = typename K::Point_3;
  using Vector_3 = typename K::Vector_3;
  using Plane_3 = typename K::Plane_3;

  VPM vpm1 = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np1, CGAL::internal_np::vertex_point), CGAL::get_property_map(CGAL::vertex_point, mesh1));
  VPM vpm2 = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np2, CGAL::internal_np::vertex_point), CGAL::get_property_map(CGAL::vertex_point, mesh2));

  using Neighbor_query1 = CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query<TriangleMesh1>;
  using Region_type1 = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region<K, TriangleMesh1>;
  using Sorting1 = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_sorting<K, TriangleMesh1, Neighbor_query1>;
  using Region_growing1 = CGAL::Shape_detection::Region_growing<Neighbor_query1, Region_type1>;

  using Neighbor_query2 = CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query<TriangleMesh2>;
  using Region_type2 = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region<K, TriangleMesh2>;
  using Sorting2 = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_sorting<K, TriangleMesh2, Neighbor_query2>;
  using Region_growing2 = CGAL::Shape_detection::Region_growing<Neighbor_query2, Region_type2>;

  static_assert(std::is_same<typename Mesh::Vertex_index, typename TriangleMesh2::Vertex_index>::value);
  static_assert(std::is_same<typename TriangleMesh1::Halfedge_index, typename TriangleMesh2::Halfedge_index>::value);
  static_assert(std::is_same<typename TriangleMesh1::Face_index, typename TriangleMesh2::Face_index>::value);
  using Vertex_index = typename TriangleMesh1::Vertex_index;
  using Halfedge_index = typename TriangleMesh1::Halfedge_index;
  using Face_index = typename TriangleMesh1::Face_index;

  std::vector<typename Region_growing1::Primitive_and_region> regions;
  std::vector<std::vector<Vertex_index> > vertices_per_region;

  std::vector<VPM> vpms = { vpm1, vpm2 };
  std::vector<Mesh*> meshes = { &mesh1, &mesh2 };

  std::vector<Plane_3> planes;

  Neighbor_query1 nq1(mesh1);
  Region_type1 rt1(mesh1, CGAL::parameters::maximum_distance(epsilon).maximum_angle(90).minimum_region_size(1));
  Sorting1 s1(mesh1, nq1);
  Region_growing1 rg1(CGAL::faces(mesh1), s1.ordered(), nq1, rt1);

  using Color = CGAL::IO::Color;
  auto plane_less = [](const Plane_3& x, const Plane_3& y) -> bool {

    if (x.a() != y.a()) return x.a() < y.a();
    if (x.b() != y.b()) return x.b() < y.b();
    if (x.c() != y.c()) return x.c() < y.c();
    return x.d() < y.d();
    };

  std::map < Plane_3, Color, decltype(plane_less) > plane_colors(plane_less);

  using Color_map = typename Mesh::template Property_map<Face_index, Color>;
  Color_map face_color1 = mesh1.template add_property_map<Face_index, Color>("color", Color(0, 0, 0)).first;

  std::vector<int> plane_index_map;
  std::vector<std::pair<std::size_t, Halfedge_index> > points; // std::size_t for mesh index

  // Reserve some memory
  plane_index_map.reserve(mesh1.number_of_halfedges() * 0.75);
  points.reserve(mesh1.number_of_halfedges() * 0.75);

  rg1.detect(boost::make_function_output_iterator(
    [&](const std::pair<Plane_3, typename Region_growing1::Region>& region) {
      // Generate a random color.
      Color color(
        static_cast<unsigned char>(rand() % 256),
        static_cast<unsigned char>(rand() % 256),
        static_cast<unsigned char>(rand() % 256));

      Plane_3 pl = region.first;

      if (CGAL::ORIGIN + pl.orthogonal_vector() < CGAL::ORIGIN + Vector_3(CGAL::NULL_VECTOR))
        pl = pl.opposite();

      color = plane_colors.emplace(pl, color).first->second;

      for (const auto& f : region.second) {
        put(face_color1, f, color);
        auto h = halfedge(f, mesh1);
        auto h2 = h;
        do {
          plane_index_map.push_back(static_cast<int>(planes.size()));
          points.push_back(std::make_pair(0, h2));
          h2 = next(h2, mesh1);
        } while (h2 != h);
      }
      planes.emplace_back(pl);
    }));

  std::size_t offset = planes.size();

  std::cout << "m1: " << planes.size() << " primitives detected" << std::endl;

  CGAL::IO::write_polygon_mesh("mesh1_segmented.ply", mesh1, CGAL::parameters::face_color_map(face_color1));

  Neighbor_query2 nq2(mesh2);
  Region_type2 rt2(mesh2, CGAL::parameters::maximum_distance(epsilon).maximum_angle(90).minimum_region_size(1));
  Sorting2 s2(mesh2, nq2);
  Region_growing2 rg2(CGAL::faces(mesh2), s2.ordered(), nq2, rt2);

  Color_map face_color2 = mesh2.template add_property_map<Face_index, Color>("color", Color(0, 0, 0)).first;

  plane_index_map.reserve(plane_index_map.size() + mesh2.number_of_halfedges() * 0.75);
  points.reserve(points.size() + mesh2.number_of_halfedges() * 0.75);

  rg2.detect(boost::make_function_output_iterator(
    [&](const std::pair<Plane_3, typename Region_growing2::Region>& region) {
      // Generate a random color.
      Color color(
        static_cast<unsigned char>(rand() % 256),
        static_cast<unsigned char>(rand() % 256),
        static_cast<unsigned char>(rand() % 256));

      Plane_3 pl = region.first;

      if (CGAL::ORIGIN + pl.orthogonal_vector() < CGAL::ORIGIN + Vector_3(CGAL::NULL_VECTOR))
        pl = pl.opposite();

      color = plane_colors.emplace(pl, color).first->second;

      for (const auto& f : region.second) {
        put(face_color2, f, color);
        auto h = halfedge(f, mesh2);
        auto h2 = h;
        do {
          plane_index_map.push_back(planes.size());
          points.push_back(std::make_pair(1, h2));
          h2 = next(h2, mesh2);
        } while (h2 != h);
      }
      planes.emplace_back(pl);
    }));

  std::cout << "m2: " << planes.size() << " primitives detected" << std::endl;

  CGAL::IO::write_polygon_mesh("mesh2_segmented.ply", mesh2, CGAL::parameters::face_color_map(face_color2));

  std::cout << plane_colors.size() << " unique planes " << unique_planes(planes).size() << std::endl;

  CGAL::Shape_regularization::Planes::regularize_planes(planes, points,
    CGAL::parameters::plane_index_map(CGAL::make_random_access_property_map(plane_index_map))
    .point_map(Point_from_halfedge_map<Point_3, std::pair<std::size_t, Halfedge_index>, Mesh, VPM>(&vpms, &meshes))
    .regularize_parallelism(true)
    .regularize_coplanarity(true)
    .maximum_angle(2)
    .maximum_offset(epsilon * 5));

  plane_colors.clear();
  for (Plane_3& pl : planes) {
    Color color(
      static_cast<unsigned char>(rand() % 256),
      static_cast<unsigned char>(rand() % 256),
      static_cast<unsigned char>(rand() % 256));

    if (CGAL::ORIGIN + pl.orthogonal_vector() < CGAL::ORIGIN + Vector_3(CGAL::NULL_VECTOR))
      pl = pl.opposite();

    color = plane_colors.emplace(pl, color).first->second;
  }

  std::cout << plane_colors.size() << " unique planes " << unique_planes(planes).size() << std::endl;

  auto region_map = rg1.region_map();
  for (const Face_index f : mesh1.faces()) {
    int r = get(region_map, f);
    if (r != -1)
      put(face_color1, f, plane_colors[planes[r]]);
  }

  CGAL::IO::write_polygon_mesh("mesh1_segmented.ply", mesh1, CGAL::parameters::face_color_map(face_color1));

  region_map = rg2.region_map();
  for (const Face_index f : mesh2.faces()) {
    int r = get(region_map, f);
    if (r != -1)
      put(face_color2, f, plane_colors[planes[r + offset]]);
  }
  CGAL::IO::write_polygon_mesh("mesh2_segmented.ply", mesh2, CGAL::parameters::face_color_map(face_color2));
}

int main(int argc, char* argv[])
{
  const std::string filename1 = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  const std::string filename2 = (argc > 2) ? argv[2] : CGAL::data_file_path("meshes/eight.off");
  const double epsilon        = (argc > 3) ? atof(argv[3]) : 0.1;

  Mesh mesh1, mesh2;
  if(!PMP::IO::read_polygon_mesh(filename1, mesh1) || !PMP::IO::read_polygon_mesh(filename2, mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  //std::vector<Mesh*> meshes = { &mesh1, &mesh2 };
  coregularize(mesh1, mesh2, epsilon);

  PMP::experimental::surface_snapping(mesh1, mesh2, epsilon);

  Mesh out;
  bool valid_union = PMP::corefine_and_compute_union(mesh1,mesh2, out);

  if(valid_union)
  {
    std::cout << "Union was successfully computed\n";
    CGAL::IO::write_polygon_mesh("snap_union.off", out, CGAL::parameters::stream_precision(17));
    return 0;
  }

  std::cout << "Union could not be computed\n";

  return 1;
}
