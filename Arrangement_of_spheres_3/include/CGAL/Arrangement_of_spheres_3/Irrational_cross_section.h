#ifndef CGAL_AOS3_IRRATIONAL_CROSS_SECTION_H
#define CGAL_AOS3_IRRATIONAL_CROSS_SECTION_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>
#include <CGAL/Arrangement_of_spheres_3/Rational_cross_section.h>

#define DPRINT(x) 
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Irrational_cross_section {
  CGAL_AOS3_TRAITS;
  typedef Irrational_cross_section CGAL_AOS3_TARG This;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  //typedef Cross_section_events CGAL_AOS3_TARG CSE;
  //typedef CGAL_AOS3_TYPENAME CS::Halfedge_handle Halfedge_handle;
  //typedef CGAL_AOS3_TYPENAME CS::Face_handle Face_handle;
public:

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

  Irrational_cross_section(const Traits &tr, CS &cs): tr_(tr),
						      cs_(cs)
						      {}

  


  /*
    If more than one face, make a vertices pass Then make another edges
    pass. For each arc, check if the location of the point differs in
    one direction. If so, compute the intersection circle of the rule
    plane of the point with the sphere. If there is no contact, then the
    edge is OK. If there is contact and the beginning/end of the circle
    surround the insertion point, then we remove the face from
    consideration.
  */



 

  /*Slice::Face_handle Slice::locate_point(const T::Sphere_point_3 & ep) {
    t_.set_temp_sphere(ep.sphere());
    return locate_point(ep, T::Key::temp_key());
    }*/



  template <class It>
  CGAL_AOS3_TYPENAME CS::Face_handle 
  locate(It b, It e,
	 CGAL_AOS3_TYPENAME Traits::Sphere_3_key k) {

    std::vector<int> locations(tr_.number_of_sphere_3s(), 0);
    std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> faces;
    std::vector<CGAL_AOS3_TYPENAME CS::Curve> edges;
    
    {
      std::cout << "Initial faces are " << std::endl;
      for (It c= b; c != e; ++c){
	cs_.write(*c,  std::cout ) << std::endl;
      }
    }

    for (It fit = b; fit != e; ++fit){
      if (!cs_.is_in_slice(*fit)) continue;
      bool ok=locate_point_check_face(k, *fit, locations/*, edges*/);
      if (ok) faces.push_back(*fit);
    }
    
    {
      std::cout << "After point check " << std::endl;
      for (unsigned int i=0; i< faces.size(); ++i){
	cs_.write(faces[i], std::cout ) << std::endl;
      }
    }

    CGAL_assertion(!faces.empty());
    if (faces.size() > 1) {
      std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> clean_faces;
      for (unsigned int i=0; i< faces.size(); ++i){
	if (locate_point_check_face_arcs(k, faces[i], locations)) {
	  clean_faces.push_back(faces[i]);
	} 
      }
      std::swap(faces, clean_faces);
    }
    
    {
      std::cout << "After arcs check " << std::endl;
      for (unsigned int i=0; i< faces.size(); ++i){
	cs_.write(faces[i], std::cout ) << std::endl;
      }
    }



    if (faces.size() > 1) {
      std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> clean_faces;
      for (unsigned int i=0; i< faces.size(); ++i){
	if (locate_point_check_face_vertices(k, faces[i])) {
	  clean_faces.push_back(faces[i]);
	} 
      }
      std::swap(faces, clean_faces);
    }
    
    {
      std::cout << "After vertices check " << std::endl;
      for (unsigned int i=0; i< faces.size(); ++i){
	cs_.write(faces[i], std::cout ) << std::endl;
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
      CGAL_assertion(0);
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
      CGAL_assertion(0);
      
    }
  
    CGAL_assertion(0);
    return CGAL_AOS3_TYPENAME CS::Face_handle();
  }

  
  CGAL_AOS3_TYPENAME CS::Face_handle 
  locate(const CGAL_AOS3_TYPENAME Traits::Sphere_3_key ep);



  CGAL_AOS3_TYPENAME CS::Halfedge_handle 
  find_rule_vertex(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &t, 
		   CGAL_AOS3_TYPENAME CS::Face_handle f,
		   CGAL_AOS3_TYPENAME CS::Curve rule) {
    CGAL_AOS3_TYPENAME CS::Halfedge_handle h;
    std::cout << "Searching for vertex for rule " << rule << " in face ";
    cs_.write(f, std::cout) << std::endl;
    try {
      h= shoot_rule(t, f, rule.key(), rule.rule_direction());
      
      //check_edge_collapse(h->prev());
    } catch (On_vertex_exception e) {
      CGAL_AOS3_TYPENAME CS::Vertex_handle v= e.vertex_handle();
      // if it is a rule in the same direction, return it, otherwise
      // pick a random edge


      if (v->point().is_rule_rule() 
	  || v->point().is_sphere_rule() 
	  && v->point().rule_coordinate()== rule.constant_coordinate()) {
	// insert on vertex
	return cs_.find_halfedge(v,f);
      } else {
	// if I am shooting up, make sure I am above the point etc.
	/*Halfedge_handle h0= sds_.find_halfedge(v,f);
	  Halfedge_handle h1= h->next();
	  if (h0->curve() == h1->curve()) {
	  bool cum=false;
	  if (rule.is_vertical()) cum = !cum;
	  if (h0->curve().arc_index() ==0 || h0->curve().arc_index() ==2) cum= !cum;
	  if (cum) {
	  h= h1;
	  } else {
	  h= h0;
	  }
	  
	  } else {
	  CGAL_assertion(0);
	  }*/
	CGAL_AOS3_TYPENAME CS::Halfedge_handle vh= v->halfedge();
	do {
	  if ( vh->curve().is_arc() && vh->face() == f) {
	    h=vh;
	    break;
	  } else if (vh->curve().is_arc() && vh->opposite()->face() == f) {
	    h= vh->opposite();
	    break;
	  }
	  vh= vh->next()->opposite();
	} while (true);
      }
    }
    //cs_.clean_edge(h);
    CGAL_AOS3_TYPENAME CS::Halfedge_handle v= cs_.insert_vertex(CGAL_AOS3_TYPENAME CS::Point(rule, 
											     h->curve()), 
								h);
    //check_edge_collapse(h);
    //CGAL_assertion(h->prev()->vertex()==v);
    CGAL_assertion(v->face() == f);
    return v;
  }




  CGAL_AOS3_TYPENAME CS::Halfedge_handle shoot_rule(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t,
						    CGAL_AOS3_TYPENAME CS::Face_handle f,
						    const CGAL_AOS3_TYPENAME Traits::Sphere_3_key & center,
						    CGAL_AOS3_TYPENAME CS::Rule_direction rd) ;


  // return comparison of point on edge of face to the shot rule on the C coordinate
  // i.e. SMALLER if the point on edge is smaller than the rule coordinate
  CGAL::Comparison_result rule_shoot_edge_vertex(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t,
						 const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt, 
						 CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						 CGAL_AOS3_TYPENAME CS::Curve hp,
						 CGAL_AOS3_TYPENAME CS::Point p,
						 CGAL_AOS3_TYPENAME CS::Curve hn) const ;

  CGAL::Comparison_result rule_shoot_compare_SS(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
						const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt, 
						CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						CGAL_AOS3_TYPENAME CS::Point point) const ;
  



  CGAL::Comparison_result rule_shoot_compare_SR(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
						const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt,
						CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						CGAL_AOS3_TYPENAME CS::Curve arc,
						CGAL_AOS3_TYPENAME Traits::Sphere_3_key orule,
						CGAL_AOS3_TYPENAME CS::Point debug_pt,
						bool arc_above) const ;


  CGAL::Comparison_result debug_rule_shoot_answer(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
						  const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt,
						  CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						  CGAL_AOS3_TYPENAME CS::Point point,
						  bool &exact) const ;

  
  
  void debug_rule_shoot_check(CGAL::Comparison_result check, 
			      CGAL::Comparison_result computed,
			      bool exact) const;


  /*bool Slice::rule_shoot_compare_if_rational_arc(const T::Sphere_point_3 & ep,
    Sds::Curve rule,
    Sds::Curve a,
    CGAL::Comparison_result &ret) const {
    CGAL::Comparison_result comp= t_.compare_sphere_centers_c(a.key(), rule.key(),
    rule.constant_coordinate());
    // NOTE could use cached locations
    if (rule.constant_coordinate()==plane_coordinate(1)){
    if (a.is_top() && comp == CGAL::LARGER) {
    ret= CGAL::LARGER;
    return true;
    } else if (!a.is_top() && comp ==CGAL::SMALLER){
    ret= CGAL::SMALLER;
    return true;
    }
    } else {
    if (a.is_right() && comp == CGAL::LARGER) {
    ret= CGAL::LARGER;
    return true;
    } else if (!a.is_right() && comp ==CGAL::SMALLER){
    ret= CGAL::SMALLER;
    return true;
    }
    }
    return false;
    }*/


  bool rule_shoot_compare_if_rational(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
				      const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt,
				      CGAL_AOS3_TYPENAME CS::Rule_direction rd,
				      CGAL_AOS3_TYPENAME CS::Point point,
				      CGAL::Comparison_result &ret) const ;





protected:

  bool locate_point_check_face(const CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
			       CGAL_AOS3_TYPENAME CS::Face_const_handle it,
			       std::vector<int> &locations) const ;





  // cache checked arcs
  bool locate_point_check_face_arcs( CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
				     CGAL_AOS3_TYPENAME CS::Face_const_handle f,
				     std::vector<int> &locations) const ;



  bool locate_point_check_face_vertices(const CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
					CGAL_AOS3_TYPENAME CS::Face_const_handle it) const ;


  int sphere_location(const CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
		      CGAL_AOS3_TYPENAME Traits::Sphere_3_key s) const ;




 




  bool behind_arc( CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
		   CGAL_AOS3_TYPENAME CS::Curve arc,
		   int location) const;


  void point_sphere_orientation(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
				CGAL_AOS3_TYPENAME Traits::Sphere_3_key sphere,
				std::vector<int> &locations
				/*,
				  std::vector<CS::Curve> &edges*/) const ;
  

  //! Take the vertex pointed to by h with the face h->face() and put a rule in it.
  CGAL_AOS3_TYPENAME CS::Halfedge_handle
  extend_rule(const CGAL_AOS3_TYPENAME Traits::Event_point_3 &t,
	      CGAL_AOS3_TYPENAME CS::Halfedge_handle h,
	      CGAL_AOS3_TYPENAME CS::Curve rule) {
    CGAL_AOS3_TYPENAME CS::Halfedge_handle v= find_rule_vertex(t,h->face(),rule);
    return cs_.split_face(rule, h, v);
  }

  /*
    walk along the rule deleting it. The last segment will be left if it hits a sphere extemum and the first segment will always be left. Other segments will necessarily go away. At the end the first segment will end in a T.
   */
  void 
  roll_back_rule(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &t,
		 CGAL_AOS3_TYPENAME CS::Halfedge_handle cur) {
    std::cout << "Rolling back, rolling back..." << cur->curve()
	      << " at " << cur->vertex()->point() << std::endl;
    std::vector<CGAL_AOS3_TYPENAME CS::Halfedge_handle> bits;
    bits.push_back(cur);
    do {
      bits.push_back(cs_.next_edge_on_rule(bits.back()));
    } while (bits.back() != CGAL_AOS3_TYPENAME CS::Halfedge_handle());
    bits.pop_back();
    
    if (bits.size() > 1) {
      
      // each internal intersection is either 3 way with only one side or 4 way
      for (unsigned int i=0; i< bits.size()-1; ++i) {
	if (cs_.degree(bits[i]->vertex())==4) {
	  // nothing to do
	} else {
	  std::cout << "Filling in vertex from edge ";
	  cs_.write(bits[i], std::cout) << " with ";
	  CGAL_AOS3_TYPENAME CS::Halfedge_handle ce= cs_.cross_edge(bits[i]); 
	  cs_.write(ce, std::cout) << std::endl;
	  CGAL_AOS3_TYPENAME CS::Face_handle f=ce->prev()->opposite()->face();
	  
	  CGAL_AOS3_TYPENAME CS::Halfedge_handle rt= extend_rule(t, 
								 cs_.find_halfedge(bits[i]->vertex(), f),
								 ce->opposite()->curve());
	  
	  
	  //std::cout << "Adding rule from ";
	  //cs_.write(find_halfedge(bits[i]->vertex(),f), std::cout) << " to ";
	  //cs_.write(rt) << std::endl;
	  //cs_.split_face(ce->opposite_curve(), find_halfedge(bits[i]->vertex(),f),
	  //      rt);				    
	}
      }
      if (bits.back()->vertex()->point().is_sphere_extremum()
	  && !bits.back()->next()->curve().is_inside()) {
	// need the back
      } else {
	CGAL_AOS3_TYPENAME CS::Vertex_handle vh= bits.back()->vertex();
	cs_.merge_faces(bits.back());
	clean_up_vertex(t, vh);
      }
      for (unsigned int i=1; i< bits.size()-1; ++i) {
	cs_.merge_faces(bits[i]);
      }
    }
 
  }							   

  void clean_up_vertex(const CGAL_AOS3_TYPENAME Traits::Event_point_3 &t,
		       CGAL_AOS3_TYPENAME CS::Vertex_handle vh) {
    if (cs_.is_redundant(vh)) {
      cs_.remove_vertex(vh->halfedge());
    }
    if (vh->point().is_sphere_extremum()) {
      CGAL_AOS3_TYPENAME Traits::Sphere_3_key k= vh->point().sphere_key();
      CGAL_AOS3_TYPENAME CS::Halfedge_handle rule, arc, c= vh->halfedge();
      do {
	if (c->curve().is_rule() && c->curve().key() == k) {
	  rule= c;
	} 
	if (!c->curve().is_rule()) {
	  arc= c;
	}
	c= c->next()->opposite();
      } while (c != vh->halfedge());
      if (rule == CGAL_AOS3_TYPENAME CS::Halfedge_handle()) {
	if (arc->curve().is_inside()) arc=arc->opposite();
	CGAL_AOS3_TYPENAME CS::Curve rule= CS::Curve::make_rule(arc->curve().key(),
								vh->point().sphere_extremum_index());
	CGAL_AOS3_TYPENAME CS::Halfedge_handle ovh
	  = find_rule_vertex(t, arc->face(),
			     rule);
	cs_.split_face(rule, cs_.find_halfedge(vh, arc->face()), ovh);
      }
    }
    
  }

  // makes rule go away
  /*
  CGAL_AOS3_TYPENAME CS::Halfedge_handle 
  rotate_rule(const CGAL_AOS3_TYPENAME Traits::Event_point_3 &ep,
	      CGAL_AOS3_TYPENAME CS::Halfedge_handle rule) {
    CGAL_precondition(rule->curve().is_rule());
    CGAL_precondition(rule->vertex()->point().is_rule_rule());
  
    std::cout << "Rotating rule " << rule->curve() << " about vertex "
	      << rule->vertex()->point() << std::endl;

    CGAL_AOS3_TYPENAME CS::Halfedge_handle oe;
    CGAL_AOS3_TYPENAME CS::Face_handle f;
    CGAL_AOS3_TYPENAME CS::Curve ec;
    if (rule->next()->curve().key() == rule->curve().key()) {
      oe= rule;
      f= rule->face();
      ec= rule->opposite()->prev()->curve();
    } else {
      oe= rule->opposite()->prev();
      f= rule->opposite()->face();
      ec= rule->next()->opposite()->curve();
    }
      
    CGAL_AOS3_TYPENAME  CS::Halfedge_handle hv= find_rule_vertex(ep, f,   ec);
    std::cout << "New edge supported by " << ec 
	      << " from " << oe->vertex()->point() << " to " 
	      << hv->vertex()->point() << std::endl;
    CGAL_AOS3_TYPENAME  CS::Halfedge_handle nh=cs_.split_face(ec, oe, hv);
    //cse_.check_edge_collapse(nh);
  
    CGAL_AOS3_TYPENAME  CS::Vertex_handle t= rule->vertex();
    std::cout << "Removing rule ";
    cs_.write(rule, std::cout);
    std::cout << std::endl;
    //cse_.clean_edge(rule);
    
    CGAL_AOS3_TYPENAME CS::Halfedge_handle th= rule->opposite()->prev();
    CGAL_assertion(th->vertex() == t);

    //cse_.check_merged_faces(rule->face(), rule->opposite()->face());
    CGAL_AOS3_TYPENAME  CS::Halfedge_handle hpa= rule->next()->opposite();
    CGAL_AOS3_TYPENAME  CS::Halfedge_handle hpb= rule->prev();
  
    cs_.merge_faces(rule);
    //cse_.check_remove_vertex(hpa);
    //cse_.check_remove_vertex(hpb);
    return hv;
    }*/



  /*CGAL_AOS3_TYPENAME CS::Face_handle  
  merge_faces(CGAL_AOS3_TYPENAME CS::Halfedge_handle h) {
    //cse_.clean_edge(h);
    //cse_.check_merged_faces(h->face(), h->opposite()->face());
    return cs_.merge_faces(h);
    }*/

  /*CGAL_AOS3_TYPENAME CS::Face_handle  
  merge_faces(CGAL_AOS3_TYPENAME CS::Vertex_handle h) {
    //cse_.clean_edge(h);
    //cse_.check_merged_faces(h->face(), h->opposite()->face());
    return cs_.merge_faces(h);
    }*/

  /*CGAL_AOS3_TYPENAME CS::Halfedge_handle 
  split_face(CGAL_AOS3_TYPENAME CS::Curve c, 
	     CGAL_AOS3_TYPENAME CS::Halfedge_handle a, 
	     CGAL_AOS3_TYPENAME CS::Halfedge_handle b) {
    //cse_.clean_edge(a);
    //cse_.clean_edge(b);
    //cse_.clean_edge(a->next());
    //cse_.clean_edge(b->next());
    CGAL_AOS3_TYPENAME CS::Halfedge_handle n= cs_.split_face(c, a, b);
    //cse_.check_edge_collapse(n);
    //cse_.check_edge_collapse(a);
    //cse_.check_edge_collapse(b);
    //cse_.check_edge_collapse(n->next());
    //cse_.check_edge_collapse(n->opposite()->next());
    return n;
    }*/

  Traits tr_;
  CS &cs_;
  //CSE //cse_;
  //CGAL_AOS3_TYPENAME Traits::Event_point_3 z_;
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include "Irrational_cross_section_impl.h"
#endif
#endif
