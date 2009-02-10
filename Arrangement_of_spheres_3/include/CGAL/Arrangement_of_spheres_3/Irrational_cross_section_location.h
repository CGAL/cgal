#ifndef CGAL_AOS3_IRRATIONAL_CROSS_SECTION_LOCATION_H
#define CGAL_AOS3_IRRATIONAL_CROSS_SECTION_LOCATION_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>
#include <CGAL/Arrangement_of_spheres_3/Rational_cross_section.h>


CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Irrational_cross_section_location {
  CGAL_AOS3_TRAITS;
  typedef Irrational_cross_section_location CGAL_AOS3_TARG This;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  //typedef Cross_section_events CGAL_AOS3_TARG CSE;
  //typedef CGAL_AOS3_TYPENAME CS::Halfedge_handle Halfedge_handle;
  //typedef CGAL_AOS3_TYPENAME CS::Face_handle Face_handle;
public:

  typedef CGAL_AOS3_TYPENAME Traits::Sphere_3_key Sphere_3_key;
  typedef CGAL_AOS3_TYPENAME Traits::Event_point_3 Event_point_3;
  typedef CGAL_AOS3_TYPENAME Traits::Sphere_point_3 Sphere_point_3;

  struct On_edge_exception {
    On_edge_exception(CGAL_AOS3_TYPENAME CS::Halfedge_handle h): h_(h){}
    CGAL_AOS3_TYPENAME CS::Halfedge_handle halfedge_handle() const {
      return h_;
    }
    CGAL_AOS3_TYPENAME CS::Halfedge_handle h_;
  };
  struct On_vertex_exception {
    On_vertex_exception(CGAL_AOS3_TYPENAME CS::Vertex_handle h): v_(h){}
    CGAL_AOS3_TYPENAME CS::Vertex_handle vertex_handle() const {
      return v_;
    }
    CGAL_AOS3_TYPENAME CS::Vertex_handle v_;
  };

  Irrational_cross_section_location(const Traits &tr, CS &cs): tr_(tr),
							       cs_(cs)
  {}

  


 

  /*Slice::Face_handle Slice::locate_point(const T::Sphere_point_3 & ep) {
    t_.set_temp_sphere(ep.sphere());
    return locate_point(ep, T::Key::temp_key());
    }*/


  //! I need the actual point for handling intersect
  template <class It, class ID, class Oit>
  void 
  incident(It b, It e,
	   ID k, Oit out) {

    std::vector<int> locations(tr_.number_of_sphere_3s(), 0);
    std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> faces;
      
    {
      CGAL_LOG(Log::LOTS, "Initial faces are " << std::endl);
      for (It c= b; c != e; ++c){
	CGAL_LOG_WRITE(Log::LOTS, cs_.write(*c,  LOG_STREAM ) << std::endl);
      }
    }

    for (It fit = b; fit != e; ++fit){
      if (!cs_.is_in_slice(*fit)) continue;
      bool ok=check_face(k, *fit, locations/*, edges*/);
      if (ok) faces.push_back(*fit);
    }
    
    {
      CGAL_LOG(Log::LOTS, "After point check " << std::endl);
      for (unsigned int i=0; i< faces.size(); ++i){
	CGAL_LOG_WRITE(Log::LOTS, cs_.write(faces[i], LOG_STREAM )
		       << std::endl);
      }
    }

    CGAL_assertion(!faces.empty());
    if (faces.size() > 1) {
      std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> clean_faces;
      for (unsigned int i=0; i< faces.size(); ++i){
	if (check_face_arcs(k, faces[i], locations)) {
	  clean_faces.push_back(faces[i]);
	} 
      }
      std::swap(faces, clean_faces);
    }
    
    {
      CGAL_LOG(Log::LOTS, "After arcs check " << std::endl);
      for (unsigned int i=0; i< faces.size(); ++i){
	CGAL_LOG_WRITE(Log::LOTS, cs_.write(faces[i], LOG_STREAM ) << std::endl);
      }
    }
    std::copy(faces.begin(), faces.end(), out);
  }


 template <class It, class ID>
  CGAL_AOS3_TYPENAME CS::Face_handle 
  locate(It b, It e,
	 ID k) {
   std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> faces;
   std::vector<CGAL_AOS3_TYPENAME CS::Curve> edges;
   
   incident(b,e,k, std::back_inserter(faces));



    if (faces.size() > 1) {
      std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> clean_faces;
      for (unsigned int i=0; i< faces.size(); ++i){
	if (check_face_vertices(k, faces[i])) {
	  clean_faces.push_back(faces[i]);
	} 
      }
      std::swap(faces, clean_faces);
    }
    
    {
      CGAL_LOG(Log::LOTS, 
	       "After vertices check " << std::endl);
      for (unsigned int i=0; i< faces.size(); ++i){
	CGAL_LOG_WRITE(Log::LOTS, cs_.write(faces[i], LOG_STREAM ) 
		       << std::endl);
      }
    }

    if (faces.size() ==1) {
      return faces[0];
    } else if (faces.size() ==2) {
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h= faces[0]->halfedge();
      do {
	if (h->opposite()->face() == faces[1]) throw On_edge_exception(h);
	h= h->next();
      } while (h != faces[0]->halfedge());
      
      cs_.write(faces[0], std::cerr) << std::endl;
      cs_.write(faces[1], std::cerr) << std::endl;
      CGAL_AOS3_TYPENAME Traits::Degeneracy_exception d;
      d.new_face(faces[0]);
      d.new_face(faces[1]);
      throw d;
    } else {
      CGAL_assertion(!faces.empty());
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h= faces[0]->halfedge();
      do {
	bool ok=true;
	for (unsigned int i=1; i< faces.size(); ++i) {
	  if (cs_.has_vertex(faces[i], h->vertex())) {
	    
	  } else { ok=false;}
	}
	if (ok)   throw On_vertex_exception(h->vertex());
	h= h->next();
      } while (h != faces[0]->halfedge());
      
      
      // must handle degeneracy for intersection
      std::cerr << "Remaining faces are " << std::endl;
      for (unsigned int i=0; i< faces.size(); ++i){
	cs_.write(faces[i], std::cerr ) << std::endl;
      }
      {
	CGAL_AOS3_TYPENAME Traits::Degeneracy_exception d;
	for (unsigned int i=0; i< faces.size(); ++i){
	  d.new_face(faces[i]);
	}
	
	throw d;
      }
    }
  
    CGAL_assertion(0);
    return CGAL_AOS3_TYPENAME CS::Face_handle();
  }





  CGAL_AOS3_TYPENAME CS::Halfedge_handle 
  find_rule_vertex(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &t, 
		   CGAL_AOS3_TYPENAME CS::Face_handle f,
		   CGAL_AOS3_TYPENAME CS::Curve rule) ;


  template <class K>
  CGAL_AOS3_TYPENAME CS::Face_handle 
  locate(const K& k){
    std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> faces;
    for (CGAL_AOS3_TYPENAME CS::Face_iterator it = cs_.faces_begin();it != cs_.faces_end(); ++it){
      faces.push_back(it);
    }
    return locate(faces.begin(), faces.end(), k);
  }


  CGAL_AOS3_TYPENAME CS::Face_handle 
  locate( CGAL_AOS3_TYPENAME Traits::Sphere_3_key k){
    std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> faces;
    for (CGAL_AOS3_TYPENAME CS::Face_iterator it = cs_.faces_begin();it != cs_.faces_end(); ++it){
      faces.push_back(it);
    }
    return locate(faces.begin(), faces.end(), 
		  CGAL_AOS3_TYPENAME Traits::Center_point_3(k, tr_.sphere_events(k).first));
  }
 

  template <class K, class Oit>
  void
  halfedges(const K &k, CGAL_AOS3_TYPENAME CS::Face_handle f, Oit out) {
    CGAL_AOS3_TYPENAME CS::Halfedge_handle c= f->halfedge();
    do {
      // too many comparisons, fix this
      Comparison_result scrx= compare_points(k, c->opposite()->vertex()->point(), 
					     plane_coordinate(0));
      Comparison_result scry= compare_points(k, c->opposite()->vertex()->point(), 
					     plane_coordinate(1));
      Comparison_result tcrx= compare_points(k, c->vertex()->point(), 
					     plane_coordinate(0));
      Comparison_result tcry= compare_points(k, c->vertex()->point(), 
					     plane_coordinate(1));
      if (c->curve().is_rule()) {
	if (c->curve().is_vertical() && (scrx == EQUAL && scry != tcry)
	    || c->curve().is_horizontal() && (scry == EQUAL && scrx != tcrx)) {
	  *out=c; ++out;
	}
      } else {
	if (scrx != scry && tcrx != tcry) {
	  Oriented_side os= tr_.oriented_side(k, c->curve().key());
	  if (os == ON_ORIENTED_BOUNDARY){
	    *out=c;
	    ++out;
	  }
	}
      }
      c=c->next();
    } while (c != f->halfedge());
  }


  template <class K>
  Comparison_result
  compare_points_c(const K &k,
		   CGAL_AOS3_TYPENAME CS::Point pt,
		   Coordinate_index C) const {
    CGAL_LOG(Log::LOTS, "Compare points called for " << k 
	     << " and " << pt << " on " << C << std::endl);
    if (pt.is_rule_rule()) {
      return tr_.compare_to_rule_c(k, pt.rule_key(C), C);
    } else if (pt.is_sphere_rule()) {
      if (pt.is_sphere_extremum()) {
	return tr_.compare_to_circle_extremum_c(k, pt.sphere_key(),
						 pt.sphere_extremum_index(),
						 C);
      } else {
	return tr_.compare_to_circle_rule_c(k, pt.sphere_key(),
					     pt.rule_key(),
					     pt.rule_constant_coordinate(),
					     pt.is_smaller(),
					     C);
      }
      
    } else {
      
      return tr_.compare_to_circle_circle_c(k, pt.sphere_key(0) , 
					      pt.sphere_key(1), 
					      C);
    }
  }

  bool
  equal_points(CGAL_AOS3_TYPENAME CS::Vertex_const_handle vh,
	       const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &pt) const {
    CGAL_LOG(Log::LOTS, "Equal points called for " << vh->point() 
	     << " and " << pt << std::endl);
    for (unsigned int i=0; i< 2; ++i) {
      CGAL_LOG(Log::LOTS, "Trying " << i << std::endl);
      Coordinate_index ci=plane_coordinate(i);
      Comparison_result cr= compare_points_c(pt, vh->point(), ci);
      if (cr != EQUAL) return false;
    }
    return true;
      
  }


protected:

  
  template <class K>
  bool check_face(const K &k,
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
	point_sphere_orientation(k, sphere, locations);
	//DPRINT(std::cout << "Testing " << h->curve() <<std::endl);
	if (!h->curve().is_compatible_location(locations[sphere.input_index()])) {
	  //DPRINT(std::cout << "Nixed by edge " << h->curve() << std::endl);
	  return false;
	} else {
	  CGAL_LOG(Log::LOTS, "Passed location for " << h->curve() 
		   << " with location " << locations[sphere.input_index()] 
		   << std::endl);
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






  template <class ID>
  bool 
  check_face_arcs(ID k,
		  CGAL_AOS3_TYPENAME CS::Face_const_handle f,
		  std::vector<int> &locations) const {
    std::set<CGAL_AOS3_TYPENAME Traits::Sphere_3_key> check_arcs;
    CGAL_AOS3_TYPENAME CS::Halfedge_const_handle h= f->halfedge();
    do {
      if (h->curve().is_arc() && !h->curve().is_inside()
	  && check_arcs.find(h->curve().key()) == check_arcs.end()){
	CGAL_LOG(Log::LOTS,  "Arc test " << k << " on " << h->curve() 
		 << std::endl);
	check_arcs.insert(h->curve().key());
	bool ba=behind_arc(k, h->curve(), 
			   locations[h->curve().key().input_index()]);
	if (ba) {
	  CGAL_LOG(Log::LOTS, "Point is behind arc " << std::endl);
	  return false;
	} else {
	  CGAL_LOG(Log::LOTS, "Point is not behind arc " << std::endl);
	}
      }
      h= h->next();
    } while (h != f->halfedge());
    return true;
  }


  template <class K>
  void point_sphere_orientation(K k,
				Sphere_3_key sphere,
				std::vector<int> &locations
				/*,
				  std::vector<CS::Curve> &edges*/) const {
    
    if (locations[sphere.input_index()]==0){
      locations[sphere.input_index()]=sphere_location(k, sphere);
      //DPRINT(std::cout << "For sphere " << sphere << " got " << SL::decode(locations[sphere]) << std::endl);
    }
  }


  template <class T>
  bool check_face_vertices(const T &t,
			   CGAL_AOS3_TYPENAME CS::Face_const_handle it) const {
    CGAL_AOS3_TYPENAME CS::Halfedge_const_handle h= it->halfedge();
    do {
      if (h->vertex()->point().is_sphere_sphere() 
	  && h->curve().key() != h->next()->curve().key()
	  && (!h->curve().is_inside() && !h->next()->curve().is_inside())) {
	// NOTE what about degeneracies?  not sure if I need to handle
	// them here
	CGAL_AOS3_TYPENAME CS::Point npt= h->vertex()->point();
	if (tr_.oriented_side_of_separating_plane(t,
						  npt.sphere_key(0),
						  npt.sphere_key(1)) == CGAL::ON_NEGATIVE_SIDE) {
	  //CGAL_assertion(0);
	  CGAL_LOG(Log::LOTS, "Face nixed by vertex " << npt << std::endl);
	  return false;
	}
      }
      h= h->next();
    } while (h != it->halfedge());
    
    return true;
  }

  /*bool check_face_vertices(const Event_point_3 & k,
    CGAL_AOS3_TYPENAME CS::Face_const_handle it) const ;*/


  template <class T>
  int sphere_location(const T &k,
		      Sphere_3_key s) const {
    CGAL::Bounded_side bs=tr_.bounded_side_of_sphere(k, s);
    
    CGAL::Comparison_result xo= tr_.compare_to_rule_c(k, s, plane_coordinate(0));
    
    CGAL::Comparison_result yo= tr_.compare_to_rule_c(k, s, plane_coordinate(1));
    CGAL_LOG(Log::LOTS, "Comparisons for " << s << " are " << bs << " " 
	     << xo << " " << yo << std::endl);
    return encode_sphere_location(bs, xo, yo);
  }


  /*int sphere_location(Sphere_3_key k,
    Sphere_3_key s) const ;*/


  static int encode_sphere_location(Bounded_side bs,
				    Comparison_result xo,
				    Comparison_result yo)  {
    int r=0;
    if (bs== CGAL::ON_BOUNDED_SIDE) {
      r |= CS::Curve::lIN_BIT;
    } else if (bs==CGAL::ON_BOUNDARY) {
      r |= CS::Curve::lIN_BIT;
      r |= CS::Curve::lOUT_BIT;
    } else {
      r |= CS::Curve::lOUT_BIT;
    }
    
    
    /*tr_.compare_sphere_center_c(s,p, 
      plane_coordinate(1));*/
    if (xo  != CGAL::SMALLER) { // LARGER
      r |= CS::Curve::lR_BIT;
    } 
    if (xo != CGAL::LARGER){ // SMALLER
      r |= CS::Curve::lL_BIT;
    }
    if (yo != CGAL::SMALLER) { // LARGER
      r |= CS::Curve::lT_BIT;
    } 
    if (yo != CGAL::LARGER) { // SMALLER
      r |= CS::Curve::lB_BIT;
    }
    return r;
  }
			     
 



  template <class T>
  bool behind_arc( const T& k,
		   CGAL_AOS3_TYPENAME CS::Curve arc,
		   int location) const {
    CGAL_precondition(!arc.is_inside());
    CGAL_precondition(arc.is_compatible_location(location));
    CGAL_AOS3_TYPENAME Traits::Coordinate_index C=arc.is_weakly_incompatible(location);
    
    if (C.is_valid()) { // we know it is opposite the center
      /*CGAL::Comparison_result cr = t_.compare_sphere_center_c(arc.key(),
	ep, C);
	if (cr == CGAL::SMALLER && arc.is_negative()
	|| cr == CGAL::LARGER && !arc.is_negative()) return false;*/
      Bounded_side bs= tr_.bounded_side_of_sphere_c(k, arc.key(), other_plane_coordinate(C));
      return bs==ON_BOUNDED_SIDE;
    } else {
      return false;
    }

  }


 
  
 
  
  Traits tr_;
  CS &cs_;
  //CSE //cse_;
  //CGAL_AOS3_TYPENAME Traits::Event_point_3 z_;
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include "Irrational_cross_section_location_impl.h"
#endif
#endif
