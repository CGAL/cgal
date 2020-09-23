#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
         <2, 3, MyTraits>::type LCC;

typedef boost::graph_traits<LCC>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<LCC>::face_descriptor     face_descriptor;
typedef boost::graph_traits<LCC>::vertex_descriptor   vertex_descriptor;
typedef boost::graph_traits<LCC>::halfedge_iterator   halfedge_iterator;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/mech-holes-shark.off";
  LCC mesh;
  CGAL::read_off(filename, mesh);

  // Incrementally fill the holes
  unsigned int nb_holes = 0;
  for ( halfedge_iterator it=halfedges(mesh).begin();
        it!=halfedges(mesh).end(); ++it)
  {
    halfedge_descriptor h=*it;
    if(is_border(h,mesh))
    {
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

  std::ofstream out("filled_LCC.off");
  out.precision(17);
  CGAL::write_off(out, mesh);

  return 0;
}
