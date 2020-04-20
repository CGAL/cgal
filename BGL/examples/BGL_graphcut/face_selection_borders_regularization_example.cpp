#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/selection.h>

#include <fstream>
#include <iostream>

using Kernel = CGAL::Simple_cartesian<double>;
using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
using Face_index = Mesh::Face_index;

int main(int argc, char** argv)
{
  std::ifstream in((argc>1) ? argv[1] : "data/blobby.off");

  if(!in)
  {
    std::cerr << "Error: could not read input file" << std::endl;
    return EXIT_FAILURE;
  }

  Mesh mesh;
  CGAL::read_off (in, mesh);

  boost::unordered_map<Face_index, bool> is_selected_map;

  // randomly select 1/3 of faces
  std::size_t nb_selected_before = 0;
  CGAL::Random rand;
  for (Face_index fi : faces(mesh))
  {
    bool selected = (rand.get_double() < 1. / 3.);
    is_selected_map[fi] = selected;
    if (selected)
      nb_selected_before ++;
  }

  std::cerr << nb_selected_before << " selected before regularization" << std::endl;

  CGAL::regularize_face_selection_borders (mesh,
                                           boost::make_assoc_property_map(is_selected_map),
                                           0.5); // using weight = 0.5

  std::size_t nb_selected_after = 0;
  for (const auto& sel : is_selected_map)
    if (sel.second)
      ++ nb_selected_after;

  std::cerr << nb_selected_after << " selected after regularization" << std::endl;

  return EXIT_SUCCESS;
}
