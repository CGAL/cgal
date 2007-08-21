#ifndef CGAL_AOS3_IRRATIONAL_CROSS_SECTION_REMOVAL_H
#define CGAL_AOS3_IRRATIONAL_CROSS_SECTION_REMOVAL_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section.h>



CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Irrational_cross_section_removal: public Irrational_cross_section CGAL_AOS3_TARG {
  CGAL_AOS3_TRAITS;
  typedef Irrational_cross_section_removal CGAL_AOS3_TARG This;
  typedef Irrational_cross_section CGAL_AOS3_TARG P;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  // typedef Cross_section_events CGAL_AOS3_TARG CSE;
  // typedef CGAL_AOS3_TYPENAME CS::Halfedge_handle Halfedge_handle;
  // typedef CGAL_AOS3_TYPENAME CS::Face_handle Face_handle;
public:


  Irrational_cross_section_removal(const Traits &tr, CS &cs): P(tr, cs){}

  
  CGAL_AOS3_TYPENAME CS::Face_handle remove_sphere(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k) {
    P::cs_.audit();
    //T::Key k= f->halfedge()->curve().key();
    
     
    // check correctness of f
    CGAL_AOS3_TYPENAME CS::Halfedge_handle rules[4];
    CGAL_AOS3_TYPENAME CS::Vertex_handle vertices[4];
 
    {
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h=P::cs_.a_halfedge(k);
      int deg=0;
      do {
	++deg;
	CGAL_assertion(h->vertex()->point().is_sphere_extremum());
	int i= h->vertex()->point().sphere_extremum_index().index();
	rules[i]= h->opposite()->prev()->opposite();
	P::roll_back_rule(P::tr_.sphere_events(k).second, rules[i]);
	vertices[i]= h->opposite()->prev()->opposite()->vertex();
	
	P::cs_.write(vertices[i], std::cout) << " is vertex " << i << std::endl;
	h= h->next();
      } while (h != P::cs_.a_halfedge(k));
      CGAL_assertion(deg==4);
    }
    // roll in each until I have a target in a face
    
   
    CGAL_AOS3_TYPENAME CS::Vertex_handle iv= rules[0]->opposite()->vertex();
    CGAL_AOS3_TYPENAME CS::Face_handle f= P::cs_.snip_out(rules, rules+4, false);
    P::cs_.delete_circle(iv);
    for (unsigned int i=0; i< 4; ++i) {
      P::clean_up_vertex(P::tr_.sphere_events(k).second, vertices[i]);
    }
   
    return f;
  }
  


  
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#if 0
#ifdef CGAL_AOS3_USE_TEMPLATES
#include "Irrational_cross_section_removal_impl.h"
#endif
#endif


#endif
