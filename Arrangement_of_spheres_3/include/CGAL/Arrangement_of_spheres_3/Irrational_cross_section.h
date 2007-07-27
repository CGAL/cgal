#ifndef CGAL_AOS3_IRRATIONAL_CROSS_SECTION_H
#define CGAL_AOS3_IRRATIONAL_CROSS_SECTION_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Irrational_cross_section {
  CGAL_AOS3_TRAITS;
  typedef Irrational_cross_section CGAL_AOS3_TARG This;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
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

  Irrational_cross_section(const Traits &tr, CS &cs, 
			   const CGAL_AOS3_TYPENAME Traits::Event_point_3 &z): tr_(tr),
									     cs_(cs),
									     z_(z){}

  


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
  locate_point(It b, It e,
	       const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & ep) {

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
      bool ok=locate_point_check_face(ep, *fit, locations/*, edges*/);
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
	if (locate_point_check_face_arcs(ep, faces[i], locations)) {
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
	if (locate_point_check_face_vertices(ep, faces[i])) {
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
  locate_point(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & ep);

private:

  bool locate_point_check_face(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &z,
			       CGAL_AOS3_TYPENAME CS::Face_const_handle it,
			       std::vector<int> &locations) const ;





  // cache checked arcs
  bool locate_point_check_face_arcs(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &ep,
				    CGAL_AOS3_TYPENAME CS::Face_const_handle f,
				    std::vector<int> &locations) const ;



  bool locate_point_check_face_vertices(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &ep,
					CGAL_AOS3_TYPENAME CS::Face_const_handle it) const ;


  int sphere_location(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3& sp,
		      CGAL_AOS3_TYPENAME Traits::Sphere_3_key s) const ;




 




  bool behind_arc(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &ep,
		  CGAL_AOS3_TYPENAME CS::Curve arc,
		  int location) const;


  void point_sphere_orientation(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &time,
				CGAL_AOS3_TYPENAME Traits::Sphere_3_key sphere,
				std::vector<int> &locations
				/*,
				  std::vector<CS::Curve> &edges*/) const ;
  



  Traits tr_;
  CS &cs_;
  CGAL_AOS3_TYPENAME Traits::Event_point_3 z_;
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include "Irrational_cross_section_impl.h"
#endif
#endif
