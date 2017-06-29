#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <boost/unordered_map.hpp>
#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <fstream>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/boost/graph/helpers.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main()
{
  Polyhedron poly;
  std::ifstream input("data/head.off");
  input >> poly;

// define my selection of faces to remove
  boost::unordered_map<Polyhedron::Face_handle, bool> is_selected_map;

  const int selection_indices[] = {501, 652, 646, 322, 328, 212, 347, 345, 352, 353, 696, 697, 698, 706, 714, 2892};
  std::set<int> index_set(&selection_indices[0], &selection_indices[0]+16);

  std::vector<Polyhedron::Face_handle> faces_to_remove;
  int index = 0;
  BOOST_FOREACH(Polyhedron::Face_handle fh, faces(poly))
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

  expand_face_selection_for_removal(poly,
                                    faces_to_remove,
                                    boost::make_assoc_property_map(is_selected_map));

  index=0;
  BOOST_FOREACH(Polyhedron::Face_handle fh, faces(poly))
  {
    if (is_selected_map[fh])
    {
      CGAL::Euler::remove_face(fh->halfedge(), poly);
      ++index;
    }
  }

  CGAL_assertion(index == 25);
  CGAL_assertion(is_valid(poly));
  return 0;
}

