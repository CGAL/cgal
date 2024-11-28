#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip> // for std::setw
#include <optional>

namespace PMP = ::CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point = K::Point_3;
using Vector = K::Vector_3;
using Plane = K::Plane_3;

using Mesh = CGAL::Surface_mesh<Point>;

using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;

bool is_reflex(const edge_descriptor e,
               const Mesh& g)
{
  auto compute_plane = [](halfedge_descriptor h, const Mesh& g) -> std::optional<Plane>
  {
    auto vpm = get(CGAL::vertex_point, g);

    std::vector<Point> face_points { get(vpm, source(h, g)),
                                     get(vpm, target(h, g)) };
    CGAL_assertion(face_points[0] != face_points[1]);

    halfedge_descriptor current = next(h, g);
    do
    {
      const Point& c = get(vpm, target(current, g));
      if(!CGAL::collinear(face_points[0], face_points[1], c))
      {
        face_points.push_back(c);
        break;
      }

      current = next(current, g);
    }
    while(current != h);

    if(face_points.size() != 3)
      return std::nullopt;

     // Create plane from 3 points
    return Plane { face_points[0], face_points[1], face_points[2] };
  };

  halfedge_descriptor h = halfedge(e, g);

  Plane plane_l, plane_r;
  if(auto o_plane_l = compute_plane(h, g)) {
    plane_l = *o_plane_l;
  } else {
    std::cerr << "Error, failed to compute left plane in is_reflex(edge)" << std::endl;
    std::exit(1);
  }

  if(auto o_plane_r = compute_plane(opposite(h, g), g)) {
    plane_r = *o_plane_r;
  } else {
    std::cerr << "Error, failed to compute right plane in is_reflex(edge)" << std::endl;
    std::exit(1);
  }

  Vector normal_l = plane_l.orthogonal_vector();
  const Point& p_src = get(CGAL::vertex_point, g, source(h, g));
  const Point& p_tgt = get(CGAL::vertex_point, g, target(h, g));
  Vector dir { p_src, p_tgt };
  CGAL_assertion(normal_l != CGAL::NULL_VECTOR);
  CGAL_assertion(dir != CGAL::NULL_VECTOR);

  Point p = p_src + CGAL::cross_product(normal_l, dir);
  return (plane_r.oriented_side(p) == CGAL::ON_POSITIVE_SIDE);
}

bool is_reflex(const vertex_descriptor v,
               const Mesh& g)
{
  if (CGAL::internal::is_isolated(v,g ))
    return false;

  for(halfedge_descriptor h : halfedges_around_target(v, g))
  {
    if(!is_reflex(edge(h, g), g)) {
        return false;
    }
  }

  return true; // all edges are reflex
}

bool is_convex(const vertex_descriptor v,
               const Mesh& g)
{
  if (CGAL::internal::is_isolated(v,g ))
    return false;

  for(halfedge_descriptor h : halfedges_around_target(v, g))
  {
    if(is_reflex(edge(h, g), g)) {
        return false;
    }
  }

  return true; // all edges are not reflex
}

std::size_t number_of_reflex_edges(const Mesh& g)
{
  std::ofstream out("/home/mrouxell/tmp/tmp.polylines.txt");
  out.precision(17);

  std::size_t cntr = 0;
  for(edge_descriptor e : edges(g))
  {
    if(is_reflex(e, g))
    {
      out << "2 " << g.point(source(e, g)) << " " << g.point(target(e, g)) << std::endl;
      ++cntr;
    }
  }

  return cntr;
}

std::size_t number_of_reflex_vertices(const Mesh& g)
{
  std::size_t cntr = 0;
  for(vertex_descriptor v : vertices(g))
    if(is_reflex(v, g))
      ++cntr;

  return cntr;
}

std::size_t number_of_non_convex_vertices(const Mesh& g)
{
  std::size_t cntr = 0;
  for(vertex_descriptor v : vertices(g))
    if(!is_convex(v, g))
      ++cntr;

  return cntr;
}

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: "
              << argv[0] << "\n"
              << "\t[output_filename] (.txt)\n"
              << std::endl;
    std::cerr << "Input format: any format readable with CGAL::IO::read_polygon_mesh()" << std::endl;
    std::cerr << "Output format: text" << std::endl;
    return EXIT_FAILURE;
  }

  const char* input_filename = argv[1];

  std::cout << "in: " << input_filename << std::endl;

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(input_filename, sm)) {
    std::cerr << "Error: failed to read input" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh: " << num_vertices(sm) << " NV " << num_faces(sm) << " NF" << std::endl;

  // number of vertices
  std::size_t nv = num_vertices(sm);

  // number of faces
  std::size_t nf = num_faces(sm);

  // connected components
  std::size_t cc = PMP::internal::number_of_connected_components(sm);

  // volume
  FT vol = PMP::volume(sm);

  // area
  FT area = PMP::area(sm);

  // reflex edges
  std::size_t reflex_edges = number_of_reflex_edges(sm);

  // reflex vertices
  std::size_t reflex_vertices = number_of_reflex_vertices(sm);

  // non-convex vertices
  std::size_t non_convex_vertices = number_of_non_convex_vertices(sm);

  // Open a file to write JSON
  std::ofstream json_file("statistics.json");
  json_file.precision(17);
  if (!json_file) {
    std::cerr << "Error: Could not open file for writing." << std::endl;
    return EXIT_FAILURE;
  }

  // Write the JSON manually
  json_file << std::string(2, ' ') << "{\n";
  json_file << std::string(4, ' ') << "\"name\": " << std::filesystem::path(input_filename).stem() << ",\n";
  json_file << std::string(4, ' ') << "\"number_of_vertices\": " << nv << ",\n";
  json_file << std::string(4, ' ') << "\"number_of_faces\": " << nf << ",\n";
  json_file << std::string(4, ' ') << "\"connected_components\": " << cc << ",\n";
  json_file << std::string(4, ' ') << "\"volume\": " << vol << ",\n";
  json_file << std::string(4, ' ') << "\"area\": " << area << ",\n";
  json_file << std::string(4, ' ') << "\"reflex_edges\": " << reflex_edges << ",\n";
  json_file << std::string(4, ' ') << "\"reflex_vertices\": " << reflex_vertices << ",\n";
  json_file << std::string(4, ' ') << "\"non_convex_vertices\": " << non_convex_vertices << "\n";
  json_file << std::string(2, ' ') << "}\n";

  return EXIT_SUCCESS;
}
