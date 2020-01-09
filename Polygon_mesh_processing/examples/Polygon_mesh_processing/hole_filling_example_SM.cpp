#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <set>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor   halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor       face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor     vertex_descriptor;

bool is_small_hole(halfedge_descriptor h, Mesh & mesh,
                   std::set<Point> & examined_points,
                   double max_hole_diam, int max_num_hole_edges)
{
  int num_hole_edges = 0;
  auto cvpm = CGAL::get_const_property_map(CGAL::vertex_point, mesh);
  CGAL::Halfedge_around_face_circulator<Mesh> circ(h, mesh), done(circ);
  CGAL::Bbox_3 hole_bbox;

  do {
    Point p = get(cvpm, target(*circ, mesh));

    if (examined_points.find(p) != examined_points.end())
      return false;
    examined_points.insert(p);

    hole_bbox += CGAL::Bbox_3(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    num_hole_edges++;

    // Exit early, to avoid unnecessary traversal of large holes
    if (num_hole_edges > max_num_hole_edges) return false;
    if (hole_bbox.xmax() - hole_bbox.xmin() > max_hole_diam) return false;
    if (hole_bbox.ymax() - hole_bbox.ymin() > max_hole_diam) return false;
    if (hole_bbox.zmax() - hole_bbox.zmin() > max_hole_diam) return false;
  } while (++circ != done);

  return true;
}

// Incrementally fill the holes that are no larger than given diameter
// and with no more than a given number of edges (if specified).

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/mech-holes-shark.off";

  // Both of these must be positive in order to be considered
  double max_hole_diam   = (argc > 2) ? boost::lexical_cast<double>(argv[2]): -1.0;
  int max_num_hole_edges = (argc > 3) ? boost::lexical_cast<int>(argv[3]) : -1;

  std::ifstream input(filename);
  Mesh mesh;
  if ( !input || !(input >> mesh) ) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  // Avoid examining a hole we studied before using a different half edge.
  std::set<Point> examined_points;

  unsigned int nb_holes = 0;
  for(halfedge_descriptor h : halfedges(mesh))
  {
    if(is_border(h,mesh))
    {

      if(max_hole_diam > 0 && max_num_hole_edges > 0 &&
         !is_small_hole(h, mesh, examined_points,
                        max_hole_diam, max_num_hole_edges))
        continue;

      std::vector<face_descriptor>  patch_facets;
      std::vector<vertex_descriptor> patch_vertices;
      bool success = std::get<0>(
        CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                  mesh,
                  h,
                  std::back_inserter(patch_facets),
                  std::back_inserter(patch_vertices),
     CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).
                  geom_traits(Kernel())) );

      std::cout << "* Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
      std::cout << "  Is fairing successful: " << success << std::endl;
      nb_holes++;
    }
  }

  std::cout << std::endl;
  std::cout << nb_holes << " holes have been filled" << std::endl;

  std::string outfile = "filled_SM.off";
  std::ofstream out(outfile.c_str());
  std::cout << "Mesh written to: " << outfile << std::endl;
  out.precision(17);
  out << mesh << std::endl;
  return 0;
}
