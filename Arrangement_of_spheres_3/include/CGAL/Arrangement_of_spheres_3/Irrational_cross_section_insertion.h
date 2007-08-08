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
  //typedef Cross_section_events CGAL_AOS3_TARG CSE;
  //typedef CGAL_AOS3_TYPENAME CS::Halfedge_handle Halfedge_handle;
  //typedef CGAL_AOS3_TYPENAME CS::Face_handle Face_handle;
public:


  Irrational_cross_section_insertion(const Traits &tr, CS &cs): P(tr, cs) {}

  
  CGAL_AOS3_TYPENAME CS::Face_handle insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
					    CGAL_AOS3_TYPENAME CS::Face_handle f) {

    // where to put the vertices to attach to
    // the halfedge points to the vertex and is on this face
    CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4];

    return finish_insert(k, f, vhs);
  }

  CGAL_AOS3_TYPENAME CS::Face_handle finish_insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
						   CGAL_AOS3_TYPENAME CS::Face_handle f,
						   CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4]) {
    std::cout << "Finishing insert into face ";
    P::cs_.write(f, std::cout);
    std::cout << "\n with vertices ";
    for (unsigned int i=0; i< 4; ++i) {
      if (vhs[i] != CGAL_AOS3_TYPENAME CS::Vertex_handle()){
	P::cs_.write(vhs[i], std::cout) << " ";
      } else {
	std::cout << "null ";
      }
    }
    std::cout << std::endl;
    //P::cs_.audit(true);
    for (unsigned int i=0; i< 4; ++i) {
      if (vhs[i] == CGAL_AOS3_TYPENAME CS::Vertex_handle()) {
	vhs[i]= find_rule_vertex(tr_.sphere_events(k).first, f,
				 CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(i)))->vertex();
      }
    }
    //P::cs_.audit(true);
    CGAL_AOS3_TYPENAME CS::Halfedge_handle rvs[4];
    for (unsigned int i=0; i< 4; ++i) {
      rvs[i]= P::cs_.find_halfedge(vhs[i], f);
    }
  
    std::cout << "Inserting target..." << std::flush;
    CGAL_AOS3_TYPENAME CS::Halfedge_handle cvs[4];
    P::cs_.new_circle(k, cvs);
    
    std::vector<CGAL_AOS3_TYPENAME CS::Halfedge_handle> out;
    const CGAL_AOS3_TYPENAME CS::Curve curves[4]={CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(0)),
						  CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(1)).other_side(),
						  CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(2)).other_side(),
						  CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(3))};
						  
    P::cs_.stitch_in(cvs, cvs+4, rvs, curves);
    std::cout << "done." << std::endl;
    for (unsigned int i=0; i< 4; ++i){
      //P::cse_.check_edge_collapse(rvs[i]->next());
      //if (P::cs_.event(rvs[i]) == CGAL_AOS3_TYPENAME CS::Event_key()) {
	//P::cse_.check_edge_collapse(rvs[i]);
      //}
      // if (P::cs_.event(rvs[i]->next()->opposite()->next()) == CGAL_AOS3_TYPENAME CS::Event_key()) {
	//P::cse_.check_edge_collapse(rvs[i]->next()->opposite()->next());
      //}
      //P::cse_.check_edge_face(rvs[i]->next()->next());
      //check_edge_collapse(rvs[i]->next()->next());
    }
    

    P::cs_.write(std::cout);
    P::cs_.audit();
    return f;    
  }

  CGAL_AOS3_TYPENAME CS::Face_handle insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
					    CGAL_AOS3_TYPENAME CS::Halfedge_handle h) {

    std::cout << "Point hit edge ";
    P::cs_.write( h, std::cout) << std::endl;
    CGAL_AOS3_TYPENAME CS::Face_handle f;
    CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4];
    if (h->curve().is_rule()) {
      P::roll_back_rule(P::tr_.sphere_events(k).first, h);
      P::roll_back_rule(P::tr_.sphere_events(k).first, h->opposite());
      int hi=h->curve().direction().index();
      int hoi=h->opposite()->curve().direction().index();
      CGAL_assertion((hi+2)%4==hoi);
      vhs[h->curve().direction().index()] = h->vertex();
      vhs[h->opposite()->curve().direction().index()]= h->opposite()->vertex();
      int hp1=(h->curve().direction().index()+1) %4;
      vhs[hp1]= find_rule_vertex(tr_.sphere_events(k).first, h->face(),
				 CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(hp1)))->vertex();
      int hop1=(h->curve().direction().index()+3) %4;
      vhs[hop1]= find_rule_vertex(tr_.sphere_events(k).first, h->opposite()->face(),
				 CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(hop1)))->vertex();
      P::cs_.audit(true);
      f= P::cs_.merge_faces(h);
    } else {
      //CGAL_AOS3_TYPENAME CS::Halfedge_handle h= e;
      if (!h->curve().is_inside()) h= h->opposite();
      f= h->face();
      std::cout << "Choosing face ";
      P::cs_.write( h, std::cout);
      std::cout << std::endl;
      int start= h->curve().arc_index();
      //P::cse_.clean_edge(h);
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
    return finish_insert(k,f,vhs);
  }

  CGAL_AOS3_TYPENAME CS::Face_handle insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
					    CGAL_AOS3_TYPENAME CS::Vertex_handle v) {
    std::cerr << "Point hit vertex " << v->point() << std::endl;
    CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4];
    CGAL_AOS3_TYPENAME CS::Face_handle f;
    if (v->point().is_rule_rule()) {
      
      std::vector<CGAL_AOS3_TYPENAME CS::Halfedge_handle> to_remove;
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h=v->halfedge()->opposite();
      CGAL_AOS3_TYPENAME CS::Face_handle missing_face;
      do {
	P::roll_back_rule(P::tr_.sphere_events(k).first, h);
	to_remove.push_back(h);
	//P::cse_.clean_edge(h);
	CGAL_AOS3_TYPENAME CS::Rule_direction rd= h->curve().direction();
	std::cout << "Curve " << h->curve() << " has direction " << rd.to_str() 
		  << " and index " << rd.index() << std::endl;
	vhs[rd.index()]= h->vertex();
	if (h->prev()->curve().direction() == h->curve().direction()) {
	  missing_face= h->face();
	}
	h= h->opposite()->next();
      } while (h != v->halfedge()->opposite()); 
      
      for (unsigned int i=0; i< 4; ++i) {
	if (vhs[i] == CGAL_AOS3_TYPENAME CS::Vertex_handle()) {
	  vhs[i]= find_rule_vertex(tr_.sphere_events(k).first, missing_face,
				   CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(i)))->vertex();

	}

      }
      P::cs_.audit(true);
      CGAL_AOS3_TYPENAME CS::Vertex_handle iv= to_remove.front()->opposite()->vertex();
      CGAL_AOS3_TYPENAME CS::Face_handle fh= P::cs_.snip_out(to_remove.begin(),
							     to_remove.end());
      P::cs_.delete_component(iv);
      /*for (unsigned int i=0; i< 4; ++i) {
	if (vhs[i] == CGAL_AOS3_TYPENAME CS::Vertex_handle()) {
	vhs[i]= find_rule_vertex(tr_.sphere_events(k).first, f,
	CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(i)))->vertex();
	}
	}*/
      return finish_insert(k, fh, vhs);
    } else if (v->point().is_sphere_rule()) {
      // move it to rule
      CGAL_AOS3_TYPENAME CS::Halfedge_handle out_rule;
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h= v->halfedge();
      do {
	if (h->curve().is_rule()){
	  return insert(k,h);
	}
	h= h->next()->opposite();
      } while (h != v->halfedge());
      CGAL_assertion(0);
      return CGAL_AOS3_TYPENAME CS::Face_handle();
    } else {
      // degeneracy
      return insert(k, v->halfedge()->face());
    }
  }
  
  
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#if 0
#ifdef CGAL_AOS3_USE_TEMPLATES
#include "Irrational_cross_section_irs_impl.h"
#endif
#endif


#endif
