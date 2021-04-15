#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO.h>

#include <boost/unordered_map.hpp>
#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <fstream>
#include <set>
#include <CGAL/boost/graph/selection.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SM;
typedef boost::graph_traits<SM>::face_descriptor face_descriptor;

int main()
{
  SM sm;
  std::ifstream input("data/head.off");
  input >> sm;

// define my selection of faces to remove
  boost::unordered_map<face_descriptor, bool> is_selected_map;

  const int selection_indices[30] = {652,18,328,698,322,212,808,353,706,869,646,352,788,696,714,796,937,2892,374,697,227,501,786,794,345,16,21,581,347,723};
  std::set<int> index_set(&selection_indices[0], &selection_indices[0]+30);

  std::vector<face_descriptor> faces_to_remove;
  int index = 0;
  for(face_descriptor fh : faces(sm))
  {
    if(index_set.count(index)==0)
      is_selected_map[fh]=false;
    else
    {
      faces_to_remove.push_back(fh);
      is_selected_map[fh]=true;
    }
    ++index;
  }

  CGAL::regularize_face_selection_borders (sm, boost::make_assoc_property_map(is_selected_map), 0.5);

  return 0;
}

