#include <CGAL/Arrangement_of_spheres_3/Slice.h>

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






bool Slice::locate_point_check_face_arcs(const T::Sphere_point_3 &ep,
					 T::Key ind,
					 Face_const_handle it,
					 std::vector<int> &locations) const {
  Halfedge_const_handle h= it->halfedge();
  do {
    if (h->curve().is_arc() && !h->curve().is_inside()){
      //std::cout << "Arc test for " << h->curve() << std::endl;
      bool ba=behind_arc(ep, ind, h->curve(), 
			 locations[h->curve().key().input_index()]);
      if (ba) {
	//std::cout << "Point is behind arc " << h->curve() << std::endl;
	return false;
      }
    }
    h= h->next();
  } while (h != it->halfedge());
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

template <class It>
Slice::Face_handle Slice::locate_point(It b, It e, const T::Sphere_point_3 & ep,
					     T::Key index) {
  /*if (CGAL::abs(ep.simple_coordinate(Coordinate_index(0))) > t_.inf()
      || CGAL::abs(ep.simple_coordinate(Coordinate_index(1))) > t_.inf()){
    std::cerr << "Coordinate out of range." << std::endl;
    CGAL_assertion(0);
    }*/
  // excessive size by 3
  std::vector<int> locations(t_.number_of_spheres(), 0);
  std::vector<Face_handle> faces;
  std::vector<Sds::Curve> edges;
  //T::Sphere_location sl= tr_.sphere_location_object(ep);
  for (It fit = b; fit != e; ++fit){
    if (!sds_.is_in_slice(fit)) continue;
    /*{
      std::cout << "Trying face ";
      write(fit, std::cout) << std::endl;
      }*/
    bool ok=locate_point_check_face(ep, fit, index, locations/*, edges*/);
    if (ok) faces.push_back(fit);
  }

  CGAL_assertion(!faces.empty());
  if (faces.size() > 1) {
    //std::cout << "simplify this " << std::endl;
    std::vector<Face_handle> clean_faces;
    for (unsigned int i=0; i< faces.size(); ++i){
      if (locate_point_check_face_arcs(ep, index, faces[i], locations)) {
	if (locate_point_check_face_vertices(ep, index, faces[i])) {
	  /*for (unsigned int i=0; i< faces.size(); ++i){
	    sds_.write_face(faces[i]->halfedge(), std::cout) << std::endl;
	  }
	  CGAL_assertion(0);*/
	  clean_faces.push_back(faces[i]);
	}

	  /*{
	    std::cout << "Face is ok ";
	    write(faces[i], std::cout) << std::endl;
	    }*/
	  //} else {
	  /*std::cout << "Face rejected on arcs ";
	    write(faces[i], std::cout) << std::endl;*/
	  //}
      } /*else {
	std::cout << "Face rejected on vertices ";
	  write(faces[i], std::cout) << std::endl;
	  }*/
    }
    std::swap(faces, clean_faces);
  }

  if (faces.size() ==1) {
    return faces[0];
  } else if (faces.size() ==2) {
    Halfedge_handle h= faces[0]->halfedge();
    do {
      if (h->opposite()->face() == faces[1]) throw On_edge_exception(h);
      h= h->next();
    } while (h != faces[0]->halfedge());

    /*if (locate_point_check_face_vertices(ep, faces[0])) {
      std::cout << "Using center line to pick face." << std::endl;
      return faces[0];
    } else {
      CGAL_assertion( locate_point_check_face_vertices(ep, faces[1]));
      std::cout << "Using center line to pick face." << std::endl;
      return faces[1];
      }*/
    //
    write(faces[0], std::cerr) << std::endl;
    write(faces[1], std::cerr) << std::endl;
    CGAL_assertion(0);
  } else {
    // what if it lands on several edges at once. Make all intersections occur before insertions
    /*
      cases:
      two edges which share a face- return that face, this is hard to detect but necessary for intersection
      search for a face which contains all common edges and return that. 

      what about a list of edges? 
    */
    CGAL_assertion(!faces.empty());
    Halfedge_handle h= faces[0]->halfedge();
    do {
      bool ok=true;
      for (unsigned int i=1; i< faces.size(); ++i) {
	if (sds_.has_vertex(faces[i], h->vertex())) {

	} else { ok=false;}
      }
      if (ok)   throw On_vertex_exception(h->vertex());
      h= h->next();
    } while (h != faces[0]->halfedge());
    for (unsigned int i=0; i< faces.size(); ++i){
      write(faces[i], std::cerr ) << std::endl;
    }
    CGAL_assertion(0);
  }
   
  CGAL_assertion(0);
  return Face_handle();
}


Slice::Face_handle 
    Slice::locate_point(const T::Sphere_point_3 & ep,
                 T::Key index) {
                   return locate_point(sds_.faces_begin(), sds_.faces_end(), ep, index);
                        }









