#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor            edge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

struct Allow_no_surface_crossing
{
  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using Point_3 = typename K::Point_3;

  CGAL::Side_of_triangle_mesh<Mesh, K> m_side_of_tmesh;

  Allow_no_surface_crossing(const Mesh& mesh)
    : m_side_of_tmesh(mesh)
  {}

  bool operator()(vertex_descriptor, Point_3 src, Point_3 tgt) const
  {
    const CGAL::Bounded_side s_src = m_side_of_tmesh(src);
    const CGAL::Bounded_side s_tgt = m_side_of_tmesh(tgt);
    return (s_src == s_tgt);
  }
};

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? std::string(argv[1]) : CGAL::data_file_path("meshes/triceratops.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return EXIT_FAILURE;
  }

  double target_edge_length = (argc > 2) ? std::stod(std::string(argv[2])) : 0.3;
  unsigned int nb_iter = (argc > 3) ? std::stoi(std::string(argv[3])) : 1;
  unsigned int nb_smoothing = (argc > 4) ? std::stoi(std::string(argv[4])) : 2;

  if (PMP::does_self_intersect(mesh))
    std::cout << "Input mesh self-intersects" << std::endl;
  else
    std::cout << "Input mesh does not self-intersect" << std::endl;

  Allow_no_surface_crossing shall_move(mesh);
  PMP::isotropic_remeshing(faces(mesh), target_edge_length, mesh,
    CGAL::parameters::number_of_iterations(nb_iter)
    .number_of_relaxation_steps(nb_smoothing)
    .allow_move_functor(shall_move)
  );

  if (PMP::does_self_intersect(mesh))
    std::cout << "Output mesh self-intersects" << std::endl;
  else
    std::cout << "Output mesh does not self-intersect" << std::endl;

  std::cout << "Remeshing done." << std::endl;
  CGAL::IO::write_polygon_mesh("out.off", mesh, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
