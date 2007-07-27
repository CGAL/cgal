#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Irrational_cross_section CGAL_AOS3_TARG ::CS::Face_handle 
Irrational_cross_section CGAL_AOS3_TARG ::locate_point(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & ep) {
    if (ep.compare(z_, sweep_coordinate()) != CGAL::EQUAL) {
      std::cerr << z_ << std::endl;
      std::cerr << ep << std::endl;
    }
    CGAL_precondition(ep.compare(z_, sweep_coordinate()) == CGAL::EQUAL);
    std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> faces;
    for (CGAL_AOS3_TYPENAME CS::Face_iterator it = cs_.faces_begin();it != cs_.faces_end(); ++it){
      faces.push_back(it);
    }
    return locate_point(faces.begin(), faces.end(), ep);
  }


CGAL_AOS3_TEMPLATE
  bool 
Irrational_cross_section CGAL_AOS3_TARG ::locate_point_check_face(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &z,
			       CGAL_AOS3_TYPENAME CS::Face_const_handle it,
			       std::vector<int> &locations) const {
    CGAL_AOS3_TYPENAME CS::Halfedge_const_handle h= it->halfedge();
    //bool finite=false;
    do {
      //if (h->curve().is_finite()) {
      // we are outside
      if (!h->curve().is_finite()) {
	if (h->curve().is_right() && h->curve().is_inside()) return false;
      } else {
	CGAL_AOS3_TYPENAME Traits::Sphere_3_key sphere= h->curve().key();
	point_sphere_orientation(z, sphere, locations);
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
CGAL_AOS3_TEMPLATE
  bool 
Irrational_cross_section CGAL_AOS3_TARG ::locate_point_check_face_arcs(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &ep,
				    CGAL_AOS3_TYPENAME CS::Face_const_handle f,
				    std::vector<int> &locations) const {
    std::set<CGAL_AOS3_TYPENAME Traits::Sphere_3_key> check_arcs;
    CGAL_AOS3_TYPENAME CS::Halfedge_const_handle h= f->halfedge();
    do {
      if (h->curve().is_arc() && !h->curve().is_inside()
	  && check_arcs.find(h->curve().key()) == check_arcs.end()){
	std::cout << "Arc test " << ep << " on " << h->curve() << std::endl;
	check_arcs.insert(h->curve().key());
	bool ba=behind_arc(ep, h->curve(), 
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


CGAL_AOS3_TEMPLATE
  bool 
Irrational_cross_section CGAL_AOS3_TARG ::locate_point_check_face_vertices(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &ep,
					CGAL_AOS3_TYPENAME CS::Face_const_handle it) const {
  
    CGAL_AOS3_TYPENAME CS::Halfedge_const_handle h= it->halfedge();
    do {
      if (h->vertex()->point().is_sphere_sphere() 
	  && h->curve().key() != h->next()->curve().key()
	  && (!h->curve().is_inside() && !h->next()->curve().is_inside())) {
	// NOTE what about degeneracies?  not sure if I need to handle
	// them here
	CGAL_AOS3_TYPENAME CS::Point npt= h->vertex()->point();
	if (tr_.oriented_side_of_separating_plane(npt.sphere_key(0),
						  npt.sphere_key(1),
						  ep) == CGAL::ON_NEGATIVE_SIDE) {
	  //CGAL_assertion(0);
	  std::cout << "Face nixed by vertex " << npt << std::endl;
	  return false;
	}
      }
      h= h->next();
    } while (h != it->halfedge());
  
    return true;
  }


CGAL_AOS3_TEMPLATE
int 
Irrational_cross_section CGAL_AOS3_TARG ::sphere_location(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3& sp,
							  CGAL_AOS3_TYPENAME Traits::Sphere_3_key s) const {
    int r=0;
    CGAL::Bounded_side bs=tr_.bounded_side_of_sphere(s,
						    sp);
    
    /*T::Event_point_3 ep=sphere_start(locate_point);
      T::Line_3 l= ep.line();
      T::Event_point_3 a(sphere(s), l);
      if (a.is_valid()) {
      T::Event_point_3 b(sphere(s), l.opposite());
      CGAL::Comparison_result ca= a.compare(ep, 2);
      CGAL::Comparison_result cb= b.compare(ep, 2);*/
    
    if (bs== CGAL::ON_BOUNDED_SIDE) {
      r |= CS::Curve::lIN_BIT;
    } else if (bs==CGAL::ON_BOUNDARY) {
      r |= CS::Curve::lIN_BIT;
      r |= CS::Curve::lOUT_BIT;
    } else {
      r |= CS::Curve::lOUT_BIT;
    }

    /*
      NOTE can make it more specific for intersection point location
    */
    
    CGAL::Comparison_result xo= tr_.compare_sphere_center_c(s, sp,
							   plane_coordinate(0));
    CGAL::Comparison_result yo= tr_.compare_sphere_center_c(s, sp, 
							    plane_coordinate(1));
    if (xo  != CGAL::LARGER) {
      r |= CS::Curve::lR_BIT;
    } 
    if (xo != CGAL::SMALLER){
      r |= CS::Curve::lL_BIT;
    }
    if (yo != CGAL::LARGER) {
      r |= CS::Curve::lT_BIT;
    } 
    if (yo != CGAL::SMALLER) {
      r |= CS::Curve::lB_BIT;
    }
    return r;
  }




 



CGAL_AOS3_TEMPLATE
bool 
Irrational_cross_section CGAL_AOS3_TARG ::behind_arc(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &ep,
		  CGAL_AOS3_TYPENAME CS::Curve arc,
		  int location) const{
    CGAL_precondition(!arc.is_inside());
    CGAL_precondition(arc.is_compatible_location(location));
    CGAL_AOS3_TYPENAME Traits::Coordinate_index C=arc.is_weakly_incompatible(location);

    if (C.is_valid()) { // we know it is opposite the center
      /*CGAL::Comparison_result cr = t_.compare_sphere_center_c(arc.key(),
	ep, C);
	if (cr == CGAL::SMALLER && arc.is_negative()
	|| cr == CGAL::LARGER && !arc.is_negative()) return false;*/
      bool b= tr_.is_over_circle_c(arc.key(), ep, C);
      return b;
    } else {
      return false;
    }
  }

CGAL_AOS3_TEMPLATE
void 
Irrational_cross_section CGAL_AOS3_TARG ::point_sphere_orientation(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &time,
								   CGAL_AOS3_TYPENAME Traits::Sphere_3_key sphere,
								   std::vector<int> &locations
				/*,
				  std::vector<CS::Curve> &edges*/) const {
  if (locations[sphere.input_index()]==0){
    locations[sphere.input_index()]=sphere_location(time, sphere);
    //DPRINT(std::cout << "For sphere " << sphere << " got " << SL::decode(locations[sphere]) << std::endl);
  }
}
  

CGAL_AOS3_END_INTERNAL_NAMESPACE
