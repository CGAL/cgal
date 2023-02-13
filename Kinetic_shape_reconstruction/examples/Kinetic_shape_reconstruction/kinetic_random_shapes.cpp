#include <list>
#include <random>
#include <vector>

#include <CGAL/Polygon_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Kinetic_shape_partition_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/point_generators_2.h>

using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;

using Kernel  = EPECK;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;
using Plane_3 = typename Kernel::Plane_3;

using Polygon_2            = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Polygon_3            = std::vector<Point_3>;

using Uniform_creator = CGAL::Creator_uniform_2<double, Point_2>;
using Point_generator = CGAL::Random_points_in_square_2<Point_2, Uniform_creator>;

using Saver = CGAL::KSR_3::Saver<Kernel>;

using IFT      = typename EPICK::FT;
using IPoint_3 = typename EPICK::Point_3;

using IPolygon_3     = std::vector<IPoint_3>;
using IPolygon_3_map = CGAL::Identity_property_map<IPolygon_3>;

using KSP = CGAL::Kinetic_shape_partition_3<EPICK>;

const std::vector<int> box_vertices_to_faces(const int i) {
  const int _vertices_to_faces[8][3] = {
    {0, 2, 4}, {1, 2, 4},
    {0, 3, 4}, {1, 3, 4},
    {0, 2, 5}, {1, 2, 5},
    {0, 3, 5}, {1, 3, 5}
  };

  const std::vector<int> faces = {
    _vertices_to_faces[i][0],
    _vertices_to_faces[i][1],
    _vertices_to_faces[i][2]
  };
  return faces;
}

const std::vector<int> box_edges_to_faces(const int i) {
  const int _faces[12][2] = {
    {0, 4}, {1, 4}, {0, 5},
    {1, 5}, {2, 4}, {3, 4},
    {2, 5}, {3, 5}, {0, 2},
    {0, 3}, {1, 2}, {1, 3}
  };

  const std::vector<int> faces = {
    _faces[i][0], _faces[i][1]
  };
  return faces;
}

const std::vector<int> box_faces_to_vertices(const int i) {
  const int _vertices[6][4] = {
    {0, 4, 6, 2}, {1, 5, 7, 3},
    {1, 5, 4, 0}, {3, 7, 6, 2},
    {1, 0, 2, 3}, {5, 4, 6, 7}
  };

  const std::vector<int> vertices = {
    _vertices[i][0], _vertices[i][1],
    _vertices[i][2], _vertices[i][3]
  };
  return vertices;
}

std::vector<int> box_faces_to_edges(const int i) {
  const int _edges[6][4] = {
    { 0, 8, 2, 9},  { 1, 10, 3, 11},
    {10, 6, 8, 4},  {11,  7, 9,  5},
    { 4, 0, 5, 1},  { 6,  2, 7,  3}
  };

  const std::vector<int> edges = {
    _edges[i][0], _edges[i][1], _edges[i][2], _edges[i][3]
  };
  return edges;
}

bool find_next_object_colliding_plane(
  const Point_3& /* pt_min */, const Point_3& /* pt_max */,
  const std::vector<Point_3>& box_corners,
  const std::vector< std::pair<std::size_t, std::size_t> >& box_edges,
  const Plane_3& plane,
  const std::vector<int>& vertices,
  const std::vector<int>& edges,
  const std::pair<bool, int>& prev_object,
  std::pair<bool, int>& next_object,
  Point_3& m) {

  for (std::size_t i = 0; i < vertices.size(); ++i) {
    const int v_i = vertices[i];
    if ((prev_object.first && prev_object.second != v_i) || (!prev_object.first)) {
      const Point_3& v = box_corners[v_i];
      if (plane.a() * v.x() + plane.b() * v.y() + plane.c() * v.z() + plane.d() == FT(0)) {
        next_object = std::make_pair(true, v_i);
        m = v;
        return true;
      }
    }
  }

  for (std::size_t i = 0; i < edges.size(); ++i) {
    const int e_i = edges[i];
    if (prev_object.first || (!prev_object.first && prev_object.second != e_i)) {
      const Point_3& s = box_corners[box_edges[e_i].first];
      const Point_3& t = box_corners[box_edges[e_i].second];

      if (
        (plane.a() * s.x() + plane.b() * s.y() + plane.c() * s.z() + plane.d()) *
        (plane.a() * t.x() + plane.b() * t.y() + plane.c() * t.z() + plane.d()) < FT(0)) {

        next_object = std::make_pair(false, e_i);

        FT x, y, z;
        if (e_i <= 3) {
          x = box_corners[box_edges[e_i].first].x(), z = box_corners[box_edges[e_i].first].z();
          y = -(plane.a() * x + plane.c() * z + plane.d()) / plane.b();
        } else if (e_i <= 7) {
          y = box_corners[box_edges[e_i].first].y(), z = box_corners[box_edges[e_i].first].z();
          x = -(plane.b() * y + plane.c() * z + plane.d()) / plane.a();
        } else {
          x = box_corners[box_edges[e_i].first].x(), y = box_corners[box_edges[e_i].first].y();
          z = -(plane.a() * x + plane.b() * y + plane.d()) / plane.c();
        }
        m = Point_3(x, y, z);
        return true;
      }
    }
  }
  return false;
}

void find_next_object_colliding_plane(
  const Point_3& pt_min, const Point_3& pt_max,
  const std::vector<Point_3>& box_corners,
  const std::vector< std::pair<std::size_t, std::size_t> >& box_edges,
  const Plane_3& plane,
  std::pair<bool, int>& next_object,
  Point_3& m) {

  std::vector<int> vertices(8, -1), edges(12, -1);
  for (std::size_t i = 0; i < vertices.size(); ++i) vertices[i] = int(i);
  for (std::size_t i = 0; i < edges.size(); ++i) edges[i] = int(i);
  std::pair<bool, int> prev_object(false, -1);
  find_next_object_colliding_plane(
    pt_min, pt_max,
    box_corners, box_edges,
    plane, vertices, edges,
    prev_object, next_object, m);
}

void construct_bounding_polygon_of_support_plane(
  const Point_3& pt_min, const Point_3& pt_max,
  const std::vector<Point_3>& box_corners,
  const std::vector< std::pair<std::size_t, std::size_t> >& box_edges,
  const Plane_3& plane,
  std::list<Point_3>& bounding_polygon,
  std::vector< std::list<int> >& bounding_faces) {

  bounding_polygon.clear();
  bounding_faces.clear();

  Point_3 m;
  std::pair<bool, int> init_object(true, -1), prev_object, curr_object;
  find_next_object_colliding_plane(
    pt_min, pt_max, box_corners, box_edges, plane, init_object, m);
  bounding_polygon.push_back(m);

  prev_object = init_object;
  int prev_face = -1;

  do {
    std::vector<int> adjacent_faces;
    if (prev_object.first) {
      adjacent_faces = box_vertices_to_faces(prev_object.second);
    } else {
      adjacent_faces = box_edges_to_faces(prev_object.second);
    }

    bool iteration_done = false;
    int curr_face;
    for (std::size_t f = 0; f < adjacent_faces.size(); ++f) {
      curr_face = adjacent_faces[f];
      if (curr_face == prev_face) continue;

      auto vertices = box_faces_to_vertices(curr_face);
      auto edges    = box_faces_to_edges(curr_face);

      iteration_done = find_next_object_colliding_plane(
        pt_min, pt_max,
        box_corners, box_edges,
        plane, vertices, edges,
        prev_object, curr_object, m);

      if (iteration_done) {
        if (curr_object != init_object) {
          bounding_polygon.push_back(m);
        }

        if (curr_object.first) {
          std::list<int> faces;
          for (std::size_t g = 0; g < adjacent_faces.size(); ++g) {
            if (adjacent_faces[g] != prev_face) {
              faces.push_back(adjacent_faces[g]);
            }
          }
          bounding_faces.push_back(faces);
        } else {
          bounding_faces.push_back(std::list<int>(1, curr_face));
        }

        prev_object = curr_object;
        prev_face   = curr_face;
        break;
      }
    }
    assert(iteration_done);
  } while (curr_object != init_object);
}

void construct_bounding_polygon_of_support_plane(
  const Point_3& pt_min, const Point_3& pt_max,
  const std::vector<Point_3>& box_corners,
  const std::vector< std::pair<std::size_t, std::size_t> >& box_edges,
  const Plane_3& plane,
  std::list<Point_3>& bounding_polygon) {

  std::vector< std::list<int> > bounding_faces;
  construct_bounding_polygon_of_support_plane(
    pt_min, pt_max,
    box_corners, box_edges,
    plane,
    bounding_polygon,
    bounding_faces);
}

void create_random_polygons(
  const std::size_t num_polygons,
  const std::size_t num_vertices,
  const double side_length,
  std::vector<Polygon_3>& polygons) {

  std::default_random_engine generator(
    std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> R(-1.0, 1.0);

  const Point_3 pt_min(-FT(1), -FT(1), -FT(1));
  const Point_3 pt_max( FT(1),  FT(1),  FT(1));
  std::vector<Point_3> box_corners;
  std::vector< std::pair<std::size_t, std::size_t> > box_edges;

  const std::size_t vertices[8][3] = {
    {0, 0, 0}, {1, 0, 0},
    {0, 1, 0}, {1, 1, 0},
    {0, 0, 1}, {1, 0, 1},
    {0, 1, 1}, {1, 1, 1}
  };
  const std::size_t edges[12][2] = {
    {0, 2}, {1, 3}, {4, 6},
    {5, 7}, {0, 1}, {2, 3},
    {4, 5}, {6, 7}, {0, 4},
    {2, 6}, {1, 5}, {3, 7}
  };

  for (std::size_t i = 0; i < 8; ++i) {
    const FT x = (vertices[i][0] == 0 ? pt_min.x() : pt_max.x());
    const FT y = (vertices[i][1] == 0 ? pt_min.y() : pt_max.y());
    const FT z = (vertices[i][2] == 0 ? pt_min.z() : pt_max.z());
    box_corners.push_back(Point_3(x, y, z));
  }

  for (std::size_t i = 0; i < 12; ++i) {
    box_edges.push_back(std::make_pair(edges[i][0], edges[i][1]));
  }

  polygons.reserve(num_polygons);
  while (polygons.size() < num_polygons) {
    const FT x_0 = static_cast<FT>(R(generator));
    const FT y_0 = static_cast<FT>(R(generator));
    const FT z_0 = static_cast<FT>(R(generator));
    const Point_3 center_ref(x_0, y_0, z_0);

    const FT a = static_cast<FT>(R(generator));
    const FT b = static_cast<FT>(R(generator));
    const FT c = static_cast<FT>(R(generator));
    const FT d = -(a * x_0 + b * y_0 + c * z_0);
    const Plane_3 plane_ref(a, b, c, d);

    std::list<Point_3>   bp_ref_3d;
    std::vector<Point_2> bp_ref_2d;

    construct_bounding_polygon_of_support_plane(
      pt_min, pt_max, box_corners, box_edges, plane_ref, bp_ref_3d);

    bp_ref_2d.reserve(bp_ref_3d.size());
    for (auto it_p = bp_ref_3d.begin(); it_p != bp_ref_3d.end(); ++it_p) {
      bp_ref_2d.push_back(plane_ref.to_2d(*it_p));
    }

    Polygon_2 bp_ref(bp_ref_2d.begin(), bp_ref_2d.end());
    if (bp_ref.orientation() == CGAL::Orientation::CLOCKWISE) {
      bp_ref.reverse_orientation();
    }

    std::vector<Point_2> ch;
    CGAL::random_convex_set_2(
      num_vertices, std::back_inserter(ch), Point_generator(side_length));

    std::vector<Point_2> k_pts;
    for (std::size_t j = 0; j < ch.size(); ++j) {
      const FT lambda = ch[j].x(), mu = ch[j].y();
      const Point_3 m = center_ref + lambda * plane_ref.base1() + mu * plane_ref.base2();
      k_pts.push_back(plane_ref.to_2d(m));
    }

    Polygon_2 k_poly(k_pts.begin(), k_pts.end());
    if (k_poly.orientation() == CGAL::Orientation::CLOCKWISE) {
      k_poly.reverse_orientation();
    }

    // std::cout << "OFF" << std::endl;
    // std::cout << "4 1 0" << std::endl;
    // for (auto it_p = k_poly.vertices_begin(); it_p != k_poly.vertices_end(); ++it_p) {
    //   std::cout << *it_p << " 0" << std::endl;
    // }
    // std::cout << "4 0 1 2 3" << std::endl;

    std::list<Polygon_with_holes_2> bp_k_intersection;
    CGAL::intersection(bp_ref, k_poly, std::back_inserter(bp_k_intersection));

    if (!bp_k_intersection.empty()) {
      for (auto it_p = bp_k_intersection.begin(); it_p != bp_k_intersection.end(); ++it_p) {
        const Polygon_2 s_poly = it_p->outer_boundary();
        std::vector<Point_3> poly_generated;
        poly_generated.reserve(s_poly.size());

        for (auto it_v = s_poly.vertices_begin(); it_v != s_poly.vertices_end(); ++it_v) {
          poly_generated.push_back(plane_ref.to_3d(*it_v));
        }
        polygons.push_back(poly_generated);
      }
    }
  }

  Saver saver;
  saver.export_polygon_soup_3(
    polygons, "rnd-polygons-" +
    std::to_string(num_polygons) + "-" + std::to_string(num_vertices));
}

int main(const int argc, const char** argv) {

  // Input.
  const std::size_t n = argc > 1 ? std::atoi(argv[1]) : 1; // number of random polygons
  const std::size_t p = argc > 2 ? std::atoi(argv[2]) : 4; // number of vertices in a polygon

  const double d = 1.5; // side of the square

  std::vector<Polygon_3> rnd_polygons;
  create_random_polygons(n, p, d, rnd_polygons);

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* input kernel: "                    << boost::typeindex::type_id<EPICK>().pretty_name()  << std::endl;
  std::cout << "* polygon kernel: "                  << boost::typeindex::type_id<Kernel>().pretty_name() << std::endl;
  std::cout << "* expected number of polygons: "     << n                                                 << std::endl;
  std::cout << "* generated number of polygons: "    << rnd_polygons.size()                               << std::endl;
  std::cout << "* number of vertices in a polygon: " << p                                                 << std::endl;
  // exit(EXIT_SUCCESS);

  IPolygon_3 input_polygon;
  std::vector<IPolygon_3> input_polygons;
  input_polygons.reserve(rnd_polygons.size());
  for (const auto& rnd_polygon : rnd_polygons) {
    input_polygon.clear();
    for (const auto& rnd_point : rnd_polygon) {
      const IFT x = static_cast<IFT>(CGAL::to_double(rnd_point.x()));
      const IFT y = static_cast<IFT>(CGAL::to_double(rnd_point.y()));
      const IFT z = static_cast<IFT>(CGAL::to_double(rnd_point.z()));
      input_polygon.push_back(IPoint_3(x, y, z));
    }
    input_polygons.push_back(input_polygon);
  }
  assert(input_polygons.size() == rnd_polygons.size());

  // Algorithm.
  KSP ksp(CGAL::parameters::verbose(true).debug(false));
  const IPolygon_3_map polygon_map;
  const unsigned int k = (argc > 3 ? std::atoi(argv[3]) : 1);
  std::cout << "* input k: " << k << std::endl;

  ksp.insert(input_polygons, polygon_map);

  ksp.initialize();

  ksp.partition(k);

  // Output.
  CGAL::Linear_cell_complex_for_combinatorial_map<3, 3> lcc;
  ksp.get_linear_cell_complex(lcc);

  std::vector<unsigned int> cells = { 0, 2, 3 }, count;
  count = lcc.count_cells(cells);

  std::cout << "For k = " << k << ":" << std::endl << " vertices: " << count[0] << std::endl << " faces: " << count[2] << std::endl << " volumes: " << count[3] << std::endl;

  std::cout << std::endl << "3D kinetic partition created in " << time << " seconds!" << std::endl << std::endl;

  return EXIT_SUCCESS;
}
