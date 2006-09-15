#include <CGAL/Arrangement_of_spheres_3/Slice.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_impl.h>
#define DPRINT(x) x












bool Slice::locate_point_check_face(const T::Sphere_point_3 &z,
				    Face_const_handle it,
				    T::Key index,
				    std::vector<int> &locations) const {
  Halfedge_const_handle h= it->halfedge();
  //bool finite=false;
  do {
    //if (h->curve().is_finite()) {
    // we are outside
    if (!h->curve().is_finite()) {
      if (h->curve().is_right() && h->curve().is_inside()) return false;
    } else {
      T::Key sphere= h->curve().key();
      point_sphere_orientation(z, index, sphere, locations);
      //DPRINT(std::cout << "Testing " << h->curve() <<std::endl);
      if (!h->curve().is_compatible_location(locations[sphere.input_index()])) {
	//DPRINT(std::cout << "Nixed by edge " << h->curve() << std::endl);
	return false;
      }
    }
    //finite=true;
    //} else {
    //std::cout << "Skipping infinite " << h->curve() << std::endl;
    //}
    h= h->next();
  } while (h != it->halfedge());
  /*if (finite) {
    std::cout << "Face " << std::endl;;
    do {
    std::cout << h->curve() << std::endl;
    h= h->next();
    } while (h != it->halfedge());
    std::cout << std::endl;
    }*/

  return true;
}





// cache checked arcs
bool Slice::locate_point_check_face_arcs(const T::Sphere_point_3 &ep,
					 T::Key ind,
					 Face_const_handle f,
					 std::vector<int> &locations) const {
  std::set<T::Key> check_arcs;
  Halfedge_const_handle h= f->halfedge();
  do {
    if (h->curve().is_arc() && !h->curve().is_inside()
	&& check_arcs.find(h->curve().key()) == check_arcs.end()){
      std::cout << "Arc test " << ep << " on " << h->curve() << std::endl;
      check_arcs.insert(h->curve().key());
      bool ba=behind_arc(ep, ind, h->curve(), 
			 locations[h->curve().key().input_index()]);
      if (ba) {
	std::cout << "Point is behind arc " << std::endl;
	return false;
      } else {
	std::cout << "Point is not behind arc " << std::endl;
      }
    }
    h= h->next();
  } while (h != f->halfedge());
  return true;
}



bool Slice::locate_point_check_face_vertices(const T::Sphere_point_3 &ep,
					     T::Key index,
					     Face_const_handle it) const {
  Halfedge_const_handle h= it->halfedge();
  do {
    if (h->vertex()->point().is_sphere_sphere() 
	&& h->curve().key() != h->next()->curve().key()
	&& (!h->curve().is_inside() && !h->next()->curve().is_inside())) {
      // NOTE what about degeneracies?  not sure if I need to handle
      // them here
      Sds::Point npt= h->vertex()->point();
      if (t_.oriented_side_of_center_plane(npt.sphere_key(0),
					   npt.sphere_key(1),
					   index) == CGAL::ON_NEGATIVE_SIDE) {
	CGAL_assertion(0);
	std::cout << "Face nixed by vertex " << npt << std::endl;
	return false;
      }
    }
    h= h->next();
  } while (h != it->halfedge());

  return true;
}





/*
  If more than one face, make a vertices pass Then make another edges
  pass. For each arc, check if the location of the point differs in
  one direction. If so, compute the intersection circle of the rule
  plane of the point with the sphere. If there is no contact, then the
  edge is OK. If there is contact and the beginning/end of the circle
  surround the insertion point, then we remove the face from
  consideration.
*/



 

Slice::Face_handle Slice::locate_point(const T::Sphere_point_3 & ep) {
  t_.set_temp_sphere(ep.sphere());
  return locate_point(ep, T::Key::temp_key());
}




Slice::Face_handle 
    Slice::locate_point(const T::Sphere_point_3 & ep,
			T::Key index) {
  return locate_point(sds_.faces_begin(), sds_.faces_end(), ep, index);
}









