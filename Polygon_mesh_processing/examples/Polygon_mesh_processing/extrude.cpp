#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <iostream>
#include <fstream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel    Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>                    Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor           vertex_descriptor;

typedef Kernel::Vector_3                                       Vector;
typedef boost::property_map<Mesh, CGAL::vertex_point_t>::type  VPMap;
typedef Mesh::template Property_map<vertex_descriptor, Vector> VNMap;

struct Bottom
{
  Bottom(VPMap pmap, VNMap nmap, double vlen)
    : pmap(pmap), nmap(nmap), vlen(vlen)
  {}


  void operator()(const vertex_descriptor& vin, const vertex_descriptor vout) const
  {
    put(pmap, vout, get(pmap, vout) - vlen* get(nmap, vin));
  }

  VPMap pmap;
  VNMap nmap;
  double vlen;
};

struct Top
{
  Top(VPMap pmap, VNMap nmap, double vlen)
    : pmap(pmap), nmap(nmap), vlen(vlen)
  {}


  void operator()(const vertex_descriptor& vin, const vertex_descriptor vout) const
  {
    put(pmap, vout, get(pmap, vout) + vlen* get(nmap, vin));
  }

  VPMap pmap;
  VNMap nmap;
  double vlen;
};

int main(int argc, char* argv[])
{
  Mesh in, out;

  std::string filename = (argc > 1) ? std::string(argv[1]) : CGAL::data_file_path("meshes/cube-ouvert.off");
  double vlen = (argc > 2) ? std::stod(argv[2]) : 0.1;

  CGAL::IO::read_polygon_mesh(filename, in);

  VNMap vnormals = in.template add_property_map<vertex_descriptor, Vector>("v:normal", CGAL::NULL_VECTOR).first;

  CGAL::Polygon_mesh_processing::compute_vertex_normals(in, vnormals);
  Bottom bottom(get(CGAL::vertex_point,out), vnormals, vlen);
  Top top(get(CGAL::vertex_point,out), vnormals, vlen);
  CGAL::Polygon_mesh_processing::extrude_mesh(in, out, bottom, top);


  filename = filename.substr(filename.find_last_of("/") + 1, filename.length() - 1);
  filename = filename.substr(0, filename.find_last_of("."));
  filename = filename + "_" + std::to_string(vlen) + ".off";
  CGAL::IO::write_polygon_mesh(filename, out, CGAL::parameters::stream_precision(17));
  return 0;
}
