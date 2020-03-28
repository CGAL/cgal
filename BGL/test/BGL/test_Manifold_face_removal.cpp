#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO.h>

#include <boost/unordered_map.hpp>
#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <fstream>
#include <set>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/boost/graph/helpers.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SM;
typedef boost::graph_traits<SM>::face_descriptor face_descriptor;

void border_cases()
{
  SM sm_ref;
  std::ifstream input("data/nm_selection_removal.off");
  input >> sm_ref;

  {
    SM sm = sm_ref;
    std::vector<SM::Face_index> faces_to_remove;
    faces_to_remove.push_back(SM::Face_index(5));
    faces_to_remove.push_back(SM::Face_index(6));
    SM::Property_map<SM::Face_index, bool> is_selected = sm.add_property_map<SM::Face_index, bool>("f:is_selected", false).first;
    is_selected[SM::Face_index(5)]=true;
    is_selected[SM::Face_index(6)]=true;
    CGAL::expand_face_selection_for_removal(faces_to_remove,
                                            sm,
                                            is_selected);
    int i=0;
    for(face_descriptor fh : sm.faces())
      if(!is_selected[fh]) ++i;
    assert(i==4);
  }

  {
    SM sm = sm_ref;
    std::vector<SM::Face_index> faces_to_remove;
    faces_to_remove.push_back(SM::Face_index(4));
    faces_to_remove.push_back(SM::Face_index(6));
    SM::Property_map<SM::Face_index, bool> is_selected = sm.add_property_map<SM::Face_index, bool>("f:is_selected", false).first;
    is_selected[SM::Face_index(4)]=true;
    is_selected[SM::Face_index(6)]=true;
    CGAL::expand_face_selection_for_removal(faces_to_remove,
                                            sm,
                                            is_selected);
    int i=0;
    for(face_descriptor fh : sm.faces())
      if(!is_selected[fh]) ++i;

    assert(i==1 || i==4); // depends on the start point
  }

  {
    SM sm = sm_ref;
    std::vector<SM::Face_index> faces_to_remove;
    faces_to_remove.push_back(SM::Face_index(4));
    faces_to_remove.push_back(SM::Face_index(5));
    faces_to_remove.push_back(SM::Face_index(6));
    SM::Property_map<SM::Face_index, bool> is_selected = sm.add_property_map<SM::Face_index, bool>("f:is_selected", false).first;
    is_selected[SM::Face_index(4)]=true;
    is_selected[SM::Face_index(5)]=true;
    is_selected[SM::Face_index(6)]=true;
    CGAL::expand_face_selection_for_removal(faces_to_remove,
                                            sm,
                                            is_selected);
    int i=0;
    for(face_descriptor fh : sm.faces())
      if(!is_selected[fh]) ++i;
    assert(i==4); // depends on the start point
  }

}

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

  std::size_t nb_input_faces = sm.number_of_faces();

  expand_face_selection_for_removal(faces_to_remove,
                                    sm,
                                    boost::make_assoc_property_map(is_selected_map));

  index=0;
  for(face_descriptor fh : faces(sm))
  {
    if (is_selected_map[fh])
    {
      CGAL::Euler::remove_face(halfedge(fh, sm), sm);
      ++index;
    }
  }

  CGAL_USE(nb_input_faces);
  assert( sm.number_of_faces()+30 < nb_input_faces);
  assert(is_valid_polygon_mesh(sm));

  border_cases();

  return 0;
}

