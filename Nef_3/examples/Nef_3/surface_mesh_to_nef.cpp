#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/OFF_to_nef_3.h>


#include <sstream>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Polygon_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron_with_indices;
typedef CGAL::Nef_polyhedron_3<K> Nef_polyhedron;

template <class PolygonMesh>
void fill_cube(PolygonMesh& pm)
{
  std::string input =
"OFF\n\
8 12 0\n\
-1 -1 -1\n\
-1 1 -1\n\
1 1 -1\n\
1 -1 -1\n\
-1 -1 1\n\
-1 1 1\n\
1 1 1\n\
1 -1 1\n\
3  0 1 3\n\
3  3 1 2\n\
3  0 4 1\n\
3  1 4 5\n\
3  3 2 7\n\
3  7 2 6\n\
3  4 0 3\n\
3  7 4 3\n\
3  6 4 7\n\
3  6 5 4\n\
3  1 5 6\n\
3  2 1 6";

  std::stringstream ss;
  ss << input;
  ss >> pm;
}

int main()
{
  // construction from a Surface_mesh
  Polygon_mesh input;
  fill_cube(input);
  Nef_polyhedron nef(input);

  Polygon_mesh output;
  std::ofstream out;

  // construction from a Polyhedron_3
  Polyhedron input_bis;
  fill_cube(input_bis);
  Nef_polyhedron nef_bis(input_bis);

  // construction from a Polyhedron_3 with indices
  Polyhedron_with_indices input_ter;
  fill_cube(input_ter);
  set_halfedgeds_items_id(input_ter);
  Nef_polyhedron nef_ter( input_ter,
                          get(CGAL::halfedge_index,input_ter),
                          get(CGAL::face_index, input_ter) );

  assert(nef==nef_bis);
  assert(nef==nef_ter);
}

