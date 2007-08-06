#ifndef CGAL_AOS3_IRRATIONAL_CROSS_SECTION_IRS_H
#define CGAL_AOS3_IRRATIONAL_CROSS_SECTION_IRS_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section.h>



CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Irrational_cross_section_insertion: public Irrational_cross_section CGAL_AOS3_TARG {
  CGAL_AOS3_TRAITS;
  typedef Irrational_cross_section_insertion CGAL_AOS3_TARG This;
  typedef Irrational_cross_section CGAL_AOS3_TARG P;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  typedef Cross_section_events CGAL_AOS3_TARG CSE;
  //typedef CGAL_AOS3_TYPENAME CS::Halfedge_handle Halfedge_handle;
  //typedef CGAL_AOS3_TYPENAME CS::Face_handle Face_handle;
public:


  Irrational_cross_section_insertion(const Traits &tr, CS &cs): P(tr, cs) {}

  
  CGAL_AOS3_TYPENAME CS::Face_handle insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
					    CGAL_AOS3_TYPENAME CS::Face_handle f) {

    // where to put the vertices to attach to
    // the halfedge points to the vertex and is on this face
    CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4];
    
    for (unsigned int i=0; i< 4; ++i) {
      if (vhs[i] == CGAL_AOS3_TYPENAME CS::Vertex_handle()) {
	vhs[i]= find_rule_vertex(tr_.sphere_events(k).first, f,
				 CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(i)))->vertex();
      }
    }
    
    return finish_insert(k, f, vhs);
  }

  CGAL_AOS3_TYPENAME CS::Face_handle finish_insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
						   CGAL_AOS3_TYPENAME CS::Face_handle f,
						   CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4]) {
    CGAL_AOS3_TYPENAME CS::Halfedge_handle rvs[4];
    for (unsigned int i=0; i< 4; ++i) {
      rvs[i]= P::cs_.find_halfedge(vhs[i], f);
    }
  
    std::cout << "Inserting target..." << std::flush;
    P::cs_.insert_target(k, rvs);
    std::cout << "done." << std::endl;
    for (unsigned int i=0; i< 4; ++i){
      P::cse_.check_edge_collapse(rvs[i]->next());
      if (P::cs_.event(rvs[i]) == CGAL_AOS3_TYPENAME CS::Event_key()) {
	P::cse_.check_edge_collapse(rvs[i]);
      }
      if (P::cs_.event(rvs[i]->next()->opposite()->next()) == CGAL_AOS3_TYPENAME CS::Event_key()) {
	P::cse_.check_edge_collapse(rvs[i]->next()->opposite()->next());
      }
      P::cse_.check_edge_face(rvs[i]->next()->next());
      //check_edge_collapse(rvs[i]->next()->next());
    }
    
    return f;    
  }

  CGAL_AOS3_TYPENAME CS::Face_handle insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
					    CGAL_AOS3_TYPENAME CS::Halfedge_handle e) {

    std::cout << "Point hit edge ";
    P::cs_.write( e, std::cout) << std::endl;
    CGAL_AOS3_TYPENAME CS::Face_handle f;
    CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4];
    if (e->curve().is_rule()) {
      f= insert_sphere_on_rule_prep(k, e, vhs);
    } else {
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h= e;
      if (!h->curve().is_inside()) h= h->opposite();
      f= h->face();
      std::cout << "Choosing face ";
      P::cs_.write( h, std::cout);
      std::cout << std::endl;
      int start= h->curve().arc_index();
      P::cse_.clean_edge(h);
      h =
	P::cs_.insert_vertex(CGAL_AOS3_TYPENAME CS::Point(CGAL_AOS3_TYPENAME CS::Curve::make_rule(k,
											       Rule_direction(start)),
						       h->curve()), h);
      vhs[start] = h->vertex();
      //h= P::cs_.find_halfedge(vhs[start], f)->next();
      h=h->next();
      vhs[(start+1)%4]=
	P::cs_.insert_vertex(CS::Point(CS::Curve::make_rule(k,
							 Rule_direction((start+1)%4)),
				    h->curve()), h)->vertex();
    }
    return finish_insert(k, f, vhs);

  }

  CGAL_AOS3_TYPENAME CS::Face_handle insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
					    CGAL_AOS3_TYPENAME CS::Vertex_handle v) {
    std::cerr << "Point hit vertex " << v->point() << std::endl;
    CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4];
    CGAL_AOS3_TYPENAME CS::Face_handle f;
    if (v->point().is_rule_rule()) {
      f= insert_sphere_on_rr_prep(k, v, vhs);
      for (unsigned int i=0; i< 4; ++i) {
	if (vhs[i] == CGAL_AOS3_TYPENAME CS::Vertex_handle()) {
	  vhs[i]= find_rule_vertex(tr_.sphere_events(k).first, f,
				   CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(i)))->vertex();
	}
      }
      return finish_insert(k, f, vhs);
    } else if (v->point().is_sphere_rule()) {
      // move it to rule
      CGAL_AOS3_TYPENAME CS::Halfedge_handle out_rule;
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h= v->halfedge();
      do {
	if (h->curve().is_rule()){
	  f= insert_sphere_on_rule_prep( k, h, vhs);
	  break;
	}
	h= h->next()->opposite();
      } while (h != v->halfedge());
      return insert(k, f);
    } else {
      // degeneracy
      CGAL_assertion(0);
    }
  }




  CGAL_AOS3_TYPENAME CS::Face_handle 
  insert_sphere_on_rule_prep(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
			     CGAL_AOS3_TYPENAME CS::Halfedge_handle h,
			     CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[]){
    std::cout << "Point hit rule " << h->curve() << std::endl;
    //Face_handle f;

    if (h->curve().is_inside()) h= h->opposite();
  
    int base=0;
    if (h->curve().is_vertical()) base=3;
    vhs[(base+1)%4] = find_rule_vertex(tr_.sphere_events(k).first, h->face(), 
				       CS::Curve::make_rule(k, 
							    Rule_direction((base+1)%4)))->vertex();
    vhs[(base+3)%4] = find_rule_vertex(tr_.sphere_events(k).first, h->opposite()->face(),
				       CS::Curve::make_rule(k,
							    Rule_direction((base+3)%4)))->vertex();
  
  
    vhs[base]= h->vertex();
    vhs[(base+2)%4]= h->opposite()->vertex();
  
    P::cse_.clean_edge(h);
    P::cs_.relabel_rule(h, CS::Curve::make_rule(k, Rule_direction(base)));
    P::cs_.relabel_rule(h->opposite(),
		     CS::Curve::make_rule(k, Rule_direction((base+2)%4)).other_side());
 

    P::cse_.clean_edge(h);
    CGAL_AOS3_TYPENAME CS::Halfedge_handle hn= h->next();
    P::cs_.merge_faces(h);

    return hn->face();
  }


  CGAL_AOS3_TYPENAME CS::Face_handle 
  insert_sphere_on_rr_prep(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
			   CGAL_AOS3_TYPENAME CS::Vertex_handle v,
			   CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[]) {
    
    {
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h=v->halfedge();
      do {
	P::cse_.clean_edge(h);
	CGAL_AOS3_TYPENAME CS::Rule_direction rd= h->opposite()->curve().rule_direction();
	vhs[rd.index()]= h->opposite()->vertex();
	h= h->opposite()->prev();
      } while (h != v->halfedge()); 
    }
    CGAL_AOS3_TYPENAME CS::Face_handle f= P::cs_.merge_faces(v);

    for (unsigned int i=0; i< 4; ++i) {
      if (vhs[i] != CGAL_AOS3_TYPENAME CS::Vertex_handle()) {
	CGAL_AOS3_TYPENAME CS::Halfedge_handle vh= P::cs_.halfedge(vhs[i], f);
	P::cse_.clean_edge(vh);
	P::cse_.clean_edge(vh->next());
      }
    }
  
    return f;
  }
  
  
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#if 0
#ifdef CGAL_AOS3_USE_TEMPLATES
#include "Irrational_cross_section_irs_impl.h"
#endif
#endif


#endif
