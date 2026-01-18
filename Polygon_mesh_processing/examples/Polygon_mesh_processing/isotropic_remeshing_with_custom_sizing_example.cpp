#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

// a sizing field that is increasing the size of edge along the y-axis
// starting at a minimum size at y-max and ending at a maximum size at
// y-min, with a linear interpolation of sizes in between the two extreme
// sizing values
struct My_sizing_field
{
  double min_size, max_size;
  double ymin, ymax;

  My_sizing_field(double min_size, double max_size, double ymin, double ymax)
    : min_size(min_size)
    , max_size(max_size)
    , ymin(ymin)
    , ymax(ymax)
  {}

  double at(K::Point_3 p) const
  {
    double y=p.y();
    return CGAL::square( (y-ymin)/(ymax-ymin) * (min_size - max_size) + max_size );
  }
  double at(const Mesh::Vertex_index v, const Mesh& mesh) const { return at(mesh.point(v)); }

  std::optional<double> is_too_long(const Mesh::Vertex_index va,
                                    const Mesh::Vertex_index vb,
                                    const Mesh& mesh) const
  {
    // TODO: no mesh as parameters?
    K::Point_3 mp = CGAL::midpoint(mesh.point(va), mesh.point(vb));
    double sql_at = at(mp);
    double sql = CGAL::squared_distance(mesh.point(va), mesh.point(vb));
    if (sql > sql_at)
      return sql / sql_at;
    return std::nullopt;
  }

  std::optional<double> is_too_short(const Mesh::Halfedge_index h,
                                     const Mesh& mesh) const
  {
    K::Point_3 mp = CGAL::midpoint(mesh.point(source(h, mesh)), mesh.point(target(h, mesh)));
    double sql_at = at(mp);
    double sql = CGAL::squared_distance(mesh.point(source(h, mesh)), mesh.point(target(h, mesh)));
    if (sql < sql_at)
      return sql / sql_at;
    return std::nullopt;
  }

  K::Point_3 split_placement(const Mesh::Halfedge_index h,
                             const Mesh& mesh) const
  {
    return CGAL::midpoint(mesh.point(source(h, mesh)), mesh.point(target(h, mesh)));
  }

  void register_split_vertex(const Mesh::Vertex_index, const Mesh&) {}
};


int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elk.off");

  Mesh mesh;
  if (!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }

  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;

  CGAL::Bbox_3 bb = PMP::bbox(mesh);
  My_sizing_field sizing_field(0.1, 30, bb.ymin(), bb.ymax());
  unsigned int nb_iter = 5;

  PMP::isotropic_remeshing(
      faces(mesh),
      sizing_field,
      mesh,
      CGAL::parameters::number_of_iterations(nb_iter)
                       .number_of_relaxation_steps(3)
      );

  CGAL::IO::write_polygon_mesh("custom_remesh_out.off", mesh, CGAL::parameters::stream_precision(17));

  std::cout << "Remeshing done." << std::endl;

  return 0;
}
