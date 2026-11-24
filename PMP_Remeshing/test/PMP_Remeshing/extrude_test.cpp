#include <CGAL/Polygon_mesh_processing/extrude.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/use.h>

#include <cassert>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SMesh;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

template<typename MAP>
struct Bot
{
  Bot(MAP map):map(map){}
  template<typename VD, typename T>
  void operator()(const T&,VD vd) const
  {
    put(map, vd, get(map, vd)+Kernel::Vector_3(-2.0,0.0,1.0));
  }
  MAP map;

};

template<typename MAP>
struct Top
{
  Top(MAP map):map(map){}

  template<typename VD, typename T>
  void operator()(const T&, VD vd) const
  {
    put(map, vd, get(map, vd)+Kernel::Vector_3(0.0,2.0,-1.0));
  }

  MAP map;
};

template <class Mesh>
std::size_t number_of_faces_of_degree(const Mesh& m, std::size_t d)
{
  std::size_t nb = 0;
  for(auto f : faces(m))
    if(CGAL::degree(f, m) == d)
      ++nb;
  return nb;
}

// `out` is the extrusion of `in`: check its combinatorics, with the faces connecting
// the bottom and the top parts being either quadrilaterals or triangles
template <class Mesh>
void check_extrusion(const Mesh& in, const Mesh& out, const bool quad_strip)
{
  std::size_t nb_border_edges = 0;
  for(auto h : halfedges(in))
    if(CGAL::is_border(h, in))
      ++nb_border_edges;

  assert(CGAL::is_valid_polygon_mesh(out));
  assert(CGAL::is_closed(out));
  assert(num_vertices(out) == 2 * num_vertices(in));

  // each border edge of `in` gives a vertical edge and a quadrilateral, or, if the strip
  // is triangulated, a vertical edge, a diagonal edge, and two triangles
  const std::size_t nb_strip_elements = (quad_strip ? 1 : 2) * nb_border_edges;
  assert(num_edges(out) == 2 * num_edges(in) + nb_strip_elements);
  assert(num_faces(out) == 2 * num_faces(in) + nb_strip_elements);

  // the bottom and the top parts keep the faces of `in` as they are
  const std::size_t nb_quads_in = number_of_faces_of_degree(in, 4);
  const std::size_t nb_triangles_in = number_of_faces_of_degree(in, 3);
  assert(number_of_faces_of_degree(out, 4) == 2 * nb_quads_in + (quad_strip ? nb_border_edges : 0));
  assert(number_of_faces_of_degree(out, 3) == 2 * nb_triangles_in + (quad_strip ? 0 : 2 * nb_border_edges));

  CGAL_USE(out, nb_strip_elements, nb_quads_in, nb_triangles_in);
}

template <class Mesh>
void test_strip(const Mesh& in)
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPMap;
  const Kernel::Vector_3 v(0.0, 0.0, -1.0);

  Mesh out;
  PMP::extrude_mesh(in, out, v);
  check_extrusion(in, out, false);

  out.clear();
  PMP::extrude_mesh(in, out, v, CGAL::parameters::default_values(),
                    CGAL::parameters::do_not_triangulate_faces(false));
  check_extrusion(in, out, false);

  out.clear();
  PMP::extrude_mesh(in, out, v, CGAL::parameters::default_values(),
                    CGAL::parameters::do_not_triangulate_faces(true));
  check_extrusion(in, out, true);

  out.clear();
  Bot<VPMap> bot(get(CGAL::vertex_point, out));
  Top<VPMap> top(get(CGAL::vertex_point, out));
  PMP::extrude_mesh(in, out, bot, top);
  check_extrusion(in, out, false);

  out.clear();
  PMP::extrude_mesh(in, out, bot, top, CGAL::parameters::default_values(),
                    CGAL::parameters::do_not_triangulate_faces(true));
  check_extrusion(in, out, true);
}

template <class Mesh>
void test_mesh(const char* filename)
{
  Mesh in, out;
  std::ifstream input(filename);

  if (!input || !(input >> in))
  {
    std::cerr << "Error: cannot read Surface Mesh : " << filename << "\n";
    assert(!CGAL::is_empty(in));
    assert(false);
    return ;
  }
  CGAL::Polygon_mesh_processing::extrude_mesh(in, out, Kernel::Vector_3(0.0, 0.0, -1.0));
  std::ofstream extruded_off("extruded.off");
  extruded_off << out;
  extruded_off.close();
  out.clear();

  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPMap;
  Bot<VPMap> bot(get(CGAL::vertex_point, out));
  Top<VPMap> top(get(CGAL::vertex_point, out));
  CGAL::Polygon_mesh_processing::extrude_mesh(in, out, bot, top);
  std::ofstream gen_extruded_off("gen_extruded.off");
  gen_extruded_off << out;
  gen_extruded_off.close();

  // triangle mesh input
  test_strip(in);

  // a single quadrilateral face
  Mesh quad;
  CGAL::make_quad(Kernel::Point_3(0,0,0), Kernel::Point_3(1,0,0),
                  Kernel::Point_3(1,1,0), Kernel::Point_3(0,1,0), quad);
  test_strip(quad);

  // a grid of quadrilaterals, with interior vertices
  Mesh grid;
  CGAL::make_grid(3, 3, grid);
  test_strip(grid);

  std::cerr << "All done." << std::endl;
}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data-remeshing/quad.off";
  test_mesh<SMesh>(filename);
  test_mesh<Polyhedron>(filename);
  return 0;
}
