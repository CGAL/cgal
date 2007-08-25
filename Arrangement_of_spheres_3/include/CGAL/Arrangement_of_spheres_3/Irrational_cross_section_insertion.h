#ifndef CGAL_AOS3_IRRATIONAL_CROSS_SECTION_IRS_H
#define CGAL_AOS3_IRRATIONAL_CROSS_SECTION_IRS_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_location.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_rules.h>



CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Irrational_cross_section_insertion {
  CGAL_AOS3_TRAITS;
  typedef Irrational_cross_section_insertion CGAL_AOS3_TARG This;
  typedef Irrational_cross_section_location CGAL_AOS3_TARG ICSL;
  typedef Irrational_cross_section_rules CGAL_AOS3_TARG ICSR;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  //typedef Cross_section_events CGAL_AOS3_TARG CSE;
  //typedef CGAL_AOS3_TYPENAME CS::Halfedge_handle Halfedge_handle;
  //typedef CGAL_AOS3_TYPENAME CS::Face_handle Face_handle;
public:


  Irrational_cross_section_insertion(const Traits &tr, CS &cs): tr_(tr), cs_(cs){}

  
  CGAL_AOS3_TYPENAME CS::Face_handle insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
					    CGAL_AOS3_TYPENAME CS::Face_handle f) ;

  CGAL_AOS3_TYPENAME CS::Face_handle finish_insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
						   CGAL_AOS3_TYPENAME CS::Face_handle f,
						   CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4]) ;

  CGAL_AOS3_TYPENAME CS::Face_handle insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
					    CGAL_AOS3_TYPENAME CS::Halfedge_handle h) ;

  CGAL_AOS3_TYPENAME CS::Face_handle insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
					    CGAL_AOS3_TYPENAME CS::Vertex_handle v) ;
  
  
private:
  Traits tr_;
  CS& cs_;
};

CGAL_AOS3_END_INTERNAL_NAMESPACE


#ifdef CGAL_AOS3_USE_TEMPLATES
#include "Irrational_cross_section_insertion_impl.h"
#endif



#endif
