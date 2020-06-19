#include <pmp/SurfaceMesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_SurfaceMesh.h>
#include <CGAL/boost/graph/properties_SurfaceMesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/boost/graph/helpers.h>

#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

typedef boost::graph_traits<pmp::SurfaceMesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<pmp::SurfaceMesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<pmp::SurfaceMesh>::face_descriptor face_descriptor;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/mech-holes-shark.off";

  pmp::SurfaceMesh sm;

  sm.read(filename);

  // Incrementally fill the holes
  unsigned int nb_holes = 0;
  for(halfedge_descriptor h : halfedges(sm))
  {
    if(CGAL::is_border(h,sm))
    {
      std::vector<face_descriptor>  patch_facets;
      std::vector<vertex_descriptor> patch_vertices;
      bool success = std::get<0>(
        CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                  sm,
                  h,
                  std::back_inserter(patch_facets),
                  std::back_inserter(patch_vertices),
                  CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, sm)).
                  geom_traits(Kernel())) );


      CGAL_assertion(CGAL::is_valid_polygon_mesh(sm));

      std::cout << "* FILL HOLE NUMBER " << ++nb_holes << std::endl;
      std::cout << "  Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
      std::cout << "  Is fairing successful: " << success << std::endl;
    }
  }
  sm.write("out.off");

  return 0;
}
