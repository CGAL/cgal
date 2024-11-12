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
#include <CGAL/Real_timer.h>

#include <iostream>
#include <string>
#include <filesystem>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

template<typename T>
struct unwrap_mesh {
  using Mesh = T;
};

template<typename T>
struct unwrap_mesh<std::reference_wrapper<T>> {
  using Mesh = T;
};

template <class Point_3, class Iterator, class MeshRange>
struct Point_from_halfedge_map {
  using key_type = Iterator;
  using value_type = Point_3;
  using reference = Point_3&;
  using category = boost::readable_property_map_tag;
  using Self = Point_from_halfedge_map<Point_3, Iterator, MeshRange>;
  using Mesh = typename unwrap_mesh<typename boost::range_value<MeshRange>::type>::Mesh;

  Point_from_halfedge_map() : meshes() {}

  Point_from_halfedge_map(MeshRange meshes) : meshes(meshes) {}

  inline friend reference
    get(Self s, key_type it) {
    Mesh& m = s.meshes[it.first];
    return m.point(target(it.second, m));
  }

  MeshRange meshes;
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

template<typename MeshRange,
  typename Intersection_kernel = CGAL::Epeck,
  typename FT = typename CGAL::Kernel_traits<typename unwrap_mesh<typename boost::range_value<MeshRange>::type>::Mesh::Point>::type::FT>
void coregularize(MeshRange meshes, FT epsilon, FT angle, FT reg_epsilon, FT reg_angle, bool save_segmented = true, bool save_regularized = true, Intersection_kernel k = Intersection_kernel()) {
  using Mesh = typename unwrap_mesh<typename boost::range_value<MeshRange>::type>::Mesh;
  using K = typename CGAL::Kernel_traits<typename Mesh::Point>::type;
  using Point_3 = typename Mesh::Point;
  using Vector_3 = typename K::Vector_3;
  using Plane_3 = typename K::Plane_3;

  using EPoint_3 = typename Intersection_kernel::Point_3;
  using ELine_3 = typename Intersection_kernel::Line_3;
  using EPlane_3 = typename Intersection_kernel::Plane_3;

  using To_exact = CGAL::Cartesian_converter<K, Intersection_kernel>;
  using From_exact = CGAL::Cartesian_converter<Intersection_kernel, K>;

  using Vertex_index = typename Mesh::Vertex_index;
  using Halfedge_index = typename Mesh::Halfedge_index;
  using Face_index = typename Mesh::Face_index;

  using Region_map = typename Mesh::template Property_map<Face_index, std::size_t>;

  using Neighbor_query = CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query<Mesh>;
  using Region_type = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region<K, Mesh>;
  using Sorting = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_sorting<K, Mesh, Neighbor_query>;
  using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type, Region_map>;

  using Color = CGAL::IO::Color;
  using Color_map = typename Mesh::template Property_map<Face_index, Color>;

  CGAL::Timer timer;
  timer.start();

  std::vector<Plane_3> planes;

  auto plane_less = [](const Plane_3& x, const Plane_3& y) -> bool {
    if (x.a() != y.a()) return x.a() < y.a();
    if (x.b() != y.b()) return x.b() < y.b();
    if (x.c() != y.c()) return x.c() < y.c();
    return x.d() < y.d();
    };
  std::map<Plane_3, Color, decltype(plane_less)> plane_colors(plane_less);

  std::vector<int> plane_index_map;
  std::vector<std::pair<std::size_t, Halfedge_index> > points; // std::size_t for mesh index
  std::vector<Region_map> region_maps(meshes.size());

  std::size_t idx = 0;
  std::vector<std::size_t> offset(meshes.size() - 1, 0);
  for (Mesh& m : meshes) {
    Neighbor_query nq(m);
    Region_type rt(m, CGAL::parameters::maximum_distance(epsilon).maximum_angle(angle).minimum_region_size(1));
    Sorting s(m, nq);

    region_maps[idx] = m.template add_property_map<Face_index, std::size_t>("region_map", std::size_t(-1)).first;

    Region_growing rg(CGAL::faces(m), s.ordered(), nq, rt, region_maps[idx]);

    Color_map face_color = m.template add_property_map<Face_index, Color>("color", Color(0, 0, 0)).first;

    plane_index_map.reserve(plane_index_map.size() + std::size_t(m.number_of_halfedges() * 0.75));
    points.reserve(points.size() + std::size_t(m.number_of_halfedges() * 0.75));

    rg.detect(boost::make_function_output_iterator(
      [&](const std::pair<Plane_3, typename Region_growing::Region>& region) {
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
          put(face_color, f, color);
          auto h = halfedge(f, m);
          auto h2 = h;
          do {
            plane_index_map.push_back(static_cast<int>(planes.size()));
            points.push_back(std::make_pair(idx, h2));
            h2 = next(h2, m);
          } while (h2 != h);
        }
        planes.emplace_back(pl);
      }));

    if (idx < offset.size())
      offset[idx] = planes.size();

    idx++;
    if (save_segmented) {
      if (idx > offset.size()) {
        std::cout << (idx - 1) << ".mesh: " << (planes.size() - offset[idx - 2]) << " primitives detected" << std::endl;
        CGAL::IO::write_polygon_mesh("mesh_" + std::to_string(idx - 1) + "_segmented.ply", m, CGAL::parameters::face_color_map(face_color).use_binary_mode(false));
      }
      else if (idx > 1) {
        std::cout << (idx - 1) << ".mesh: " << (offset[idx - 1] - offset[idx - 2]) << " primitives detected" << std::endl;
        CGAL::IO::write_polygon_mesh("mesh_" + std::to_string(idx - 1) + "_segmented.ply", m, CGAL::parameters::face_color_map(face_color).use_binary_mode(false));
      }
      else {
        std::cout << "0.mesh: " << offset[0] << " primitives detected" << std::endl;
        CGAL::IO::write_polygon_mesh("mesh_0_segmented.ply", m, CGAL::parameters::face_color_map(face_color).use_binary_mode(false));
      }
    }
  }

  std::cout << plane_colors.size() << " unique planes" << std::endl;

  FT t1 = timer.time();
  std::cout << t1 << "s" <<  std::endl;

  CGAL::Shape_regularization::Planes::regularize_planes(planes, points,
    CGAL::parameters::plane_index_map(CGAL::make_random_access_property_map(plane_index_map))
    .point_map(Point_from_halfedge_map<Point_3, std::pair<std::size_t, Halfedge_index>, MeshRange>(meshes))
    .regularize_parallelism(true)
    .regularize_coplanarity(true)
    .maximum_angle(reg_angle)
    .maximum_offset(reg_epsilon));

  std::cout << (timer.time() - t1) << "s" << std::endl;

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

  std::cout << plane_colors.size() << " unique planes" << std::endl;

  // Projection of vertices
  typedef CGAL::Halfedge_around_target_circulator<Mesh> Halfedge_around_target_circulator;
  idx = 0;
  From_exact from_exact;
  To_exact to_exact;
  for (Mesh& m : meshes) {
    for (Vertex_index v : m.vertices()) {
      std::set<Plane_3, decltype(plane_less)> pl(plane_less);
      Halfedge_around_target_circulator hatc(halfedge(v, m), m), done(hatc);
      do {
        Face_index f = face(*hatc, m);
        std::size_t r = get(region_maps[idx], f);

        if (r != std::size_t(-1)) {
          if (idx == 0)
            CGAL_assertion(r < offset[0]);
          else
            CGAL_assertion(r < (planes.size() - offset[idx - 1]));
        }

        if (r != std::size_t(-1)) {
          if (idx > 0)
            r += offset[idx - 1];

          pl.insert(planes[r]);
        }

        hatc++;
      } while (hatc != done);

      if (!pl.empty()) {
        if (pl.size() == 1)
          m.point(v) = (*pl.begin()).projection(m.point(v));
        else if (pl.size() == 2) {
          auto it = pl.begin();
          const EPlane_3& first = to_exact(*it++);
          const EPlane_3& second = to_exact(*it++);
          const auto inter = CGAL::intersection(first, second);
          CGAL_assertion(bool(inter));

          // An empty intersection can only happen at sharp features and a too strong regularization.
          if (!inter)
            continue;

          ELine_3 l;
          if (CGAL::assign(l, inter)) {
            m.point(v) = from_exact(l.projection(to_exact(m.point(v))));
          }
          else { // Should not happen as planes should be unique
            CGAL_assertion(false);
            EPlane_3 plane;
            if (CGAL::assign(plane, inter))
              m.point(v) = from_exact(plane.projection(to_exact(m.point(v))));
          }
        }
        else if (pl.size() == 3) {
          auto it = pl.begin();
          const EPlane_3& first = to_exact(*it++);
          const EPlane_3& second = to_exact(*it++);
          const EPlane_3& third = to_exact(*it++);

          const auto inter = CGAL::intersection(first, second, third);
          EPoint_3 p;
          ELine_3 l;

          if (CGAL::assign(p, inter))
            m.point(v) = from_exact(p);
          else if (CGAL::assign(l, inter))
            m.point(v) = from_exact(l.projection(to_exact(m.point(v))));
        }
        else {
          std::cout << "v" << std::size_t(v) << " more than 3 incident planes!" << std::endl;
        }
      }
    }
    idx++;
  }

  if (save_regularized) {
    idx = 0;
    for (Mesh& m : meshes) {
      Color_map face_color = m.template add_property_map<Face_index, Color>("color", Color(0, 0, 0)).first;
      std::size_t shift = (idx == 0) ? 0 : offset[idx - 1];
      for (const Face_index f : m.faces()) {
        std::size_t r = get(region_maps[idx], f);
        if (r != std::size_t(-1))
          put(face_color, f, plane_colors[planes[r + shift]]);
      }

      CGAL::IO::write_polygon_mesh("mesh_" + std::to_string(idx) + "_regularized.ply", m, CGAL::parameters::face_color_map(face_color).use_binary_mode(false));
      idx++;
    }
  }
}

std::vector<Mesh> load_folder(std::string folder) {
  std::vector<Mesh> meshes;
  for (const auto& entry : std::filesystem::directory_iterator(folder)) {
    std::cout << entry.path() << std::endl;
    Mesh mesh;
    if (PMP::IO::read_polygon_mesh(entry.path().string(), mesh))
      meshes.emplace_back(mesh);
  }

  return meshes;
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

  std::vector<std::reference_wrapper<Mesh>> meshes = { std::ref(mesh1), std::ref(mesh2) };
  //std::vector<Mesh> meshes = load_folder("C:/Data/Unite");
  coregularize(meshes, epsilon, 90.0, epsilon, 0.1, true, true, CGAL::Epick());

  // PMP::experimental::surface_snapping(mesh1, mesh2, epsilon);

  std::ofstream("reg_m1.off") << std::setprecision(17) << mesh1;
  std::ofstream("reg_m2.off") << std::setprecision(17) << mesh2;

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
