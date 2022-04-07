#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#include <iostream>
#include <string>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point;

typedef CGAL::Surface_mesh<Point>                           Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;



int main()
{
// create a mesh with many connected connected components that should
// be reoriented to define a valid closed mesh
  const std::string filename = CGAL::data_file_path("meshes/elephant.off");

  Mesh mesh;
  CGAL::IO::read_polygon_mesh(filename, mesh);

  // turn the mesh into a triangle soup, duplicating all the vertices and shuffling orientations
  std::vector<Point> points;
  std::vector< std::array<std::size_t, 3> > triangles;
  triangles.reserve(faces(mesh).size());
  points.reserve(3*triangles.size());
  for (Mesh::Face_index f : mesh.faces())
  {
    Mesh::Halfedge_index h = mesh.halfedge(f);
    std::size_t s = points.size();
    points.push_back(mesh.point(source(h,mesh)));
    points.push_back(mesh.point(target(h,mesh)));
    points.push_back(mesh.point(target(mesh.next(h),mesh)));
    triangles.push_back( {s, s+1, s+2} );
    std::random_shuffle(triangles.back().begin(), triangles.back().end());
  }

  // load the soup into the mesh;
  mesh.clear();
  PMP::polygon_soup_to_polygon_mesh(points, triangles, mesh);

// now check how face orientation should be reversed
  auto fbm = mesh.add_property_map<Mesh::Face_index, bool>("fbm", false).first;
  PMP::connected_components_compatible_orientations(mesh, fbm);

// now reverse orientation of faces with bit 1
  std::vector<Mesh::Face_index> faces_to_reverse;
  for (Mesh::Face_index f : mesh.faces())
    if (get(fbm, f))
      faces_to_reverse.push_back(f);
  PMP::reverse_face_orientations(faces_to_reverse, mesh);

// now stitch borders to close the mesh
  PMP::stitch_borders(mesh);

  assert(CGAL::is_closed(mesh));
  std::ofstream("out.off") << mesh;

  return 0;
}
