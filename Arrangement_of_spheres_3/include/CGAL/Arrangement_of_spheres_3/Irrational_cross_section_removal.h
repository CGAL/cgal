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
  //typedef Cross_section_events CGAL_AOS3_TARG CSE;
  //typedef CGAL_AOS3_TYPENAME CS::Halfedge_handle Halfedge_handle;
  //typedef CGAL_AOS3_TYPENAME CS::Face_handle Face_handle;
public:


  Irrational_cross_section_removal(const Traits &tr, CS &cs): P(tr, cs){}

  
  CGAL_AOS3_TYPENAME CS::Face_handle remove_sphere(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k) {
    P::cs_.audit();
    //T::Key k= f->halfedge()->curve().key();
    
    CGAL_AOS3_TYPENAME CS::Face_handle f;
    
    f= P::cs_.a_halfedge(k)->face();
    
    // check correctness of f
    CGAL_AOS3_TYPENAME CS::Halfedge_handle rules[4];
    CGAL_AOS3_TYPENAME CS::Halfedge_handle h=f->halfedge();
    int deg=0;
    do {
      ++deg;
      CGAL_assertion(h->curve().key()==k);
      CGAL_assertion(h->curve().is_arc());
      CGAL_assertion(h->curve().is_inside());
      /*Halfedge_handle rule= h->opposite()->next();
	int index= rule->curve().rule_index();*/
      /*CGAL_assertion(rules[index]
	== Halfedge_handle());
	rules[index] =rule;*/
      h= h->next();
    } while (h != f->halfedge());
    CGAL_assertion(deg==4);
    // roll in each until I have a target in a face
  
    CGAL_AOS3_TYPENAME Traits::Sphere_3_key keys[4];
    //Halfedge_handle vs[4];
    for (unsigned int i=0; i< 4; ++i){
      rules[i]= P::cs_.rule_halfedge(k,Rule_direction(i))->opposite();
      //keys[i]= P::roll_back_rule(tr_.sphere_events(k).second, rules[i]);
      /*if (keys[i].is_valid()) {
      vs[i]= rules[i]->opposite()->prev();
      CGAL_assertion(vs[i]->vertex() == rules[i]->vertex());
      }*/
      //P::cse_.clean_edge(rules[i]);
    }
    
    //Qt_examiner_viewer_2 *qt= new Qt_examiner_viewer_2();
    //*qt_debug_examiner_viewer_2__ << Layer(0);
    //draw_rz(qt_debug_examiner_viewer_2__, CGAL::to_double(sim_->current_time()) + .1);
    //qt_debug_examiner_viewer_2__->show_everything();
    //qt_debug_examiner_viewer_2__->show();
    //*qt_debug_examiner_viewer_2__ << std::flush;

    for (unsigned int i=0; i< 4; ++i){
      CGAL_AOS3_TYPENAME CS::Halfedge_handle hi= P::cs_.rule_halfedge(k, 
								   CGAL_AOS3_TYPENAME CS::Rule_direction(i))->next();
      /*if (sds_.event(hi) != Event_key()) {
      std::cerr << "ERROR " << sds_.event(hi) 
		<< " on edge ";
      sds_.write(hi, std::cerr);
      std::cerr << std::endl;
      }*/
      //P::cse_.clean_edge(hi->opposite()->prev());
    }
    
    //P::cse_.check_merged_faces(P::cs_.rule_halfedge(k, Rule_direction(0))->face(),
    //P::cs_.rule_halfedge(k, Rule_direction(2))->face());
    //P::cse_.check_merged_faces(P::cs_.rule_halfedge(k, Rule_direction(1))->face(),
    //		    P::cs_.rule_halfedge(k, Rule_direction(3))->face());
    //return Face_handle();
    CGAL_AOS3_TYPENAME CS::Halfedge_handle vertices[4];
    

    // HERE
    //CGAL_AOS3_TYPENAME CS::Face_handle fr= P::cs_.remove_target(rules, vertices);
  
    // I guess I don't actually have to extend the rules any more
    // it doesn't  make anything invalid and makes the structure simpler. 
    // Unless there is only halfedge between the two circles
    bool has_split=false;
    for (unsigned int i=0; i< 4; ++i){
      if (keys[i].is_valid()) {
	Rule_direction rd((i+2)%4);
	CGAL_AOS3_TYPENAME CS::Curve c=CGAL_AOS3_TYPENAME CS::Curve::make_rule(keys[i], rd);
	//CGAL_assertion(vertices[i]->vertex()->point().is_sphere_rule());
	CGAL_AOS3_TYPENAME CS::Halfedge_handle hv;
	if (!has_split) {
	  CGAL_assertion(vertices[rd.index()] != CGAL_AOS3_TYPENAME CS::Halfedge_handle());
	  std::cout << "using orphaned extremum " 
		    << vertices[i]->vertex()->point() 
		    << " " << vertices[rd.index()]->vertex()->point() << std::endl;
	  // I need to handle this
	  hv= vertices[rd.index()];
	  //vertices[(i+2)%4]=Halfedge_handle();
	  keys[rd.index()]= CGAL_AOS3_TYPENAME Traits::Sphere_3_key();
	  has_split=true;
	} else {
	  std::cout << "Fixing orphaned extremum " 
		    << vertices[i]->vertex()->point() << std::endl;
	  // I need to handle this
	  // I could just insert it in the new edge since I know it goes there
	  hv= P::find_rule_vertex(P::tr_.sphere_events(k).second, vertices[i]->face(),c);
	  //vertices[i]->vertex()->point().replace_rule(c);
	}
	CGAL_AOS3_TYPENAME CS::Halfedge_handle nh;
	if (i==0 || i == 3) {
	  nh=cs_.split_face(c, hv,vertices[i]);
	} else {
	  nh=cs_.split_face(c, vertices[i], hv);
	}
	//cse_.check_edge_collapse(nh);
      } 
    }
    for (unsigned int i=0; i< 4; ++i){ 
      /*if (P::cse_.check_remove_vertex(vertices[i]) != vertices[i]) {
	if (P::cs_.is_redundant(h->vertex())
	    || P::cs_.is_redundant(h->opposite()->vertex())) {
	  // skipping since it is adjacent to another vertex which will be removed
	  // seems like there should be a better way
	  std::cout << "Skipping ephemeral edge ";
	  P::cs_.write(h, std::cout);
	  std::cout << std::endl;
	} else {
	  //P::cse_.check_edge_collapse(h);
	}
	}*/
    }

    //P::cse_.check_reduced_face(fr);

    //audit();
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
