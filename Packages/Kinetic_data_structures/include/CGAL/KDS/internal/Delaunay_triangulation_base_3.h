#ifndef KINETIC_DELAUNAY_BASE_3_H
#define KINETIC_DELAUNAY_BASE_3_H

#include <CGAL/KDS/basic.h>

#include <CGAL/KDS/internal/Triangulation_helper_3.h>
#include <CGAL/KDS/internal/triangulation_helpers_3.h>
//#include <CGAL/KDS/internal/Delaunay_cache_3.h>

// STL
#include <map>
#include <iterator>
#include <ostream>
#include <iostream>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

template <class KD, class Root_stack>
class Delaunay_event_base_3 {
public:
  Delaunay_event_base_3(const Root_stack &s,
			KD *kdel):s_(s),
				  kdel_(kdel){
  }
  //! Default constructor
  /*!  This really should not be used, but I need it for the
    Simulator::event() call, due to the apparent gcc compiler bug.
  */
  Delaunay_event_base_3(){

  }
  const Root_stack& root_stack() const {
    return s_;
  }
protected:
  const Root_stack s_;
  KD *kdel_;
};


template <class KD, class Root_stack, class Edge>
class Delaunay_3_edge_flip_event: public Delaunay_event_base_3<KD, Root_stack>  {
public:
  typedef Delaunay_event_base_3<KD, Root_stack>   P;
  Delaunay_3_edge_flip_event(const Root_stack &s,
			     const Edge &e,
			     KD *kdel):Delaunay_event_base_3<KD, Root_stack>(s, kdel), e_(e){
#ifndef NDEBUG
    o_= edge_point(e_,0);
    d_= edge_point(e_,1);
#endif
  }
  void process(const typename Root_stack::Root&){
    P::kdel_->flip(e_);
  }
  
  static typename KD::Point_key edge_point(const Edge &e, int i) {
    return vertex_of_edge(e, i)->point();
  }

  void write(std::ostream &out) const {
    out << "Flip ";
    internal::write_edge(e_, out);
#ifndef NDEBUG
    out << "(" << o_ << d_<<")" << std::flush;
    CGAL_postcondition(o_== vertex_of_edge(e_,static_cast<int>(0))->point());
    CGAL_postcondition(d_== vertex_of_edge(e_,1)->point());
#endif
  }
  const Root_stack& root_stack() const {
    return P::s_;
  }
protected:
  const Edge e_;
#ifndef NDEBUG
  typename KD::Point_key o_, d_;
#endif
};

template <class B, class C, class D>
std::ostream& operator<<(std::ostream &out, const Delaunay_3_edge_flip_event<B, C, D> &e){
  e.write(out);
  return out;
}

template <class KD, class Root_stack, class Facet>
class Delaunay_3_facet_flip_event:  public Delaunay_event_base_3<KD, Root_stack>   {
public:
  typedef Delaunay_event_base_3<KD, Root_stack>   P;
  Delaunay_3_facet_flip_event(const Root_stack &s, 
			      const Facet &e,
			      KD *kdel):  Delaunay_event_base_3<KD, Root_stack>(s, kdel),e_(e){
#ifndef NDEBUG
    a_= vertex_of_facet(e_,0)->point();
    b_= vertex_of_facet(e_,1)->point();
    c_= vertex_of_facet(e_,2)->point();
#endif
  }
  void process(const typename Root_stack::Root&){
    P::kdel_->flip(e_);
  }
  void write(std::ostream &out) const {
    out << "Flip ";
    write_facet(e_, out);
#ifndef NDEBUG
    out << "(" << a_ << b_<<c_ << ")";
    CGAL_postcondition(a_== vertex_of_facet(e_,0)->point());
    CGAL_postcondition(b_== vertex_of_facet(e_,1)->point());
    CGAL_postcondition(c_== vertex_of_facet(e_,2)->point());
#endif
  }
protected:
  const Facet e_;
#ifndef NDEBUG
  typename KD::Point_key a_, b_, c_;
#endif
};

template <class T, class C, class D>
std::ostream& operator<<(std::ostream &out, const Delaunay_3_facet_flip_event<T, C, D> &e){
  e.write(out);
  return out;
}


template <class TraitsT, class Visitor> 
class Delaunay_triangulation_base_3;

template <class TraitsT, class Visitor> 
std::ostream &operator<<(std::ostream &out, 
			 const Delaunay_triangulation_base_3<TraitsT, Visitor> &dt);

// Common base class for Delaunay and regular triangulations
template <class TraitsT, class Visitor> 
class Delaunay_triangulation_base_3  {
public:
  typedef Delaunay_triangulation_base_3<TraitsT, Visitor> This;

  // KDS typedefs
  typedef typename TraitsT::Simulator Simulator;
  typedef typename TraitsT::Moving_point_table Moving_object_table;
  typedef typename TraitsT::Kinetic_kernel Kinetic_kernel;
  //typedef typename Simulator::Time Time;
  typedef typename Moving_object_table::Key Point_key;
  typedef typename Moving_object_table::Data Point;
  typedef typename Simulator::Event_key Event_key;
  typedef typename Simulator::Root_stack Root_stack;
  //typedef typename Simulator::NT NT;

  // Delaunay typedefs
  typedef Triangulation_helper_3<typename TraitsT::Triangulation> Tri;

  // Triangulation members
  typedef typename Tri::Facet Facet;
  typedef typename Tri::Edge Edge;
  typedef typename Tri::Facet_circulator Facet_circulator;
  typedef typename Tri::Cell_circulator Cell_circulator;
  typedef typename Tri::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Tri::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tri::Edge_iterator Edge_iterator;
  //typedef typename Tri::Facet_iterator Facet_iterator;  
  typedef typename Tri::Vertex_handle Vertex_handle;
  typedef typename Tri::Cell_handle Cell_handle;
  typedef typename Tri::All_cells_iterator All_cells_iterator;
  typedef typename Tri::All_edges_iterator All_edges_iterator;
  typedef typename Tri::All_facets_iterator All_facets_iterator;
  typedef typename Tri::Finite_vertices_iterator Finite_vertices_iterator;

  // accessory types
  typedef typename TraitsT::Facet_flip Facet_flip_event;
  typedef typename Facet_flip_event::P Event_base;
  typedef typename TraitsT::Edge_flip Edge_flip_event;
  typedef Tri Triangulation;

  Delaunay_triangulation_base_3(TraitsT tr, Visitor v): tr_(tr),
							triangulation_(tr.instantaneous_kernel_object()),
							soc_(tr.side_of_oriented_sphere_3_object()),
							o3_(tr.orientation_3_object()),
							v_(v){
    if (0) print();
    has_certificates_=false;
  }
 

  const Visitor& visitor() const {
    return v_;
  }
  Visitor& visitor() {
    return v_;
  }

  //! Just write the objects in order;
  template <class Stream>
  void write(Stream &out) const {
    //out << triangulation_;
    for (All_cells_iterator cit= triangulation_.all_cells_begin(); cit != triangulation_.all_cells_end();
	 ++cit){
      Cell_handle h= cit;
      for (unsigned int i=0; i<4; ++i){
	out << h->vertex(i)->point() << " ";
      }
      out << "--";
      for (unsigned int i=0; i<4; ++i){
	Facet f(h, i);
	out << triangulation_.label(f) << " ";
      }
      out << std::endl;
    }
    for (All_edges_iterator eit = triangulation_.all_edges_begin(); eit != triangulation_.all_edges_end(); ++eit){
      triangulation_.write_labeled_edge(*eit, out);
      if (triangulation_.has_degree_3(*eit)){
	//out << " " << triangulation_.label(*eit);
      } else {
	if (has_event(*eit)){
	  out << " ?? " ;//<< triangulation_.label(*eit);
	}
	//CGAL_assertion(triangulation_.label(*eit)== Event_key::null());
      }
      out << std::endl;
    }
    for (All_facets_iterator eit = triangulation_.all_facets_begin(); eit != triangulation_.all_facets_end(); ++eit){
      triangulation_.write_labeled_facet(*eit, out);
      if (!triangulation_.has_degree_3(*eit)){
	//out << " " << triangulation_.label(*eit);
      } else {
	//CGAL_assertion(triangulation_.label(*eit)== Event_key::null());
	if (has_event(*eit)){
	  out << " ?? ";//<< triangulation_.label(*eit);
	}
      }
      out << std::endl;
    }
  }






  bool print() const {
    write(std::cout);
    return true;
  }



  /*!
    Delete all incident face and edge events (have to push edge events on to facets).
    Insert point.
    Set edges to never fails. Set outside facets.
  */
  Vertex_handle push_vertex(Point_key k, Cell_handle c){
    for (unsigned int i=0; i< 4; ++i){
      Facet f(c, i);
      if (has_event(f)){
	simulator()->delete_event(triangulation_.label(f));
	triangulation_.set_label(f, Event_key());
      }
      for (unsigned int j=0; j<i; ++j){
	Edge e(c, i, j);
	if (has_event(e)){
	  make_not_edge_flip(e, c);
	}
      }
    }
    handle_delete_cell(c);
    v_.delete_cells(&c, &c+1);
    // into max dim simplex?
    Vertex_handle vh=triangulation_.tds().insert_in_cell( c);
    vh->set_point(k);
    set_vertex_handle(k, vh);
    std::vector<Cell_handle> ics;
    triangulation_.incident_cells(vh, std::back_insert_iterator<std::vector<Cell_handle> >(ics));
    for (unsigned int j=0; j< ics.size(); ++j){
      Cell_handle cc= ics[j];
      handle_new_cell(cc);
    }
    v_.new_cells(ics.begin(), ics.end());
    for (unsigned int j=0; j< ics.size(); ++j){
      Cell_handle cc= ics[j];
      for (int i=0; i< 3; ++i){
	Edge e= triangulation_.edge_around_vertex(cc, vh, i);
	if (internal::edge_label(e)== typename Simulator::Event_key()){
	  internal::set_edge_label(triangulation(), e, simulator()->null_event());
	}
      }
      for (int i=0; i< 4; ++i){
	Facet f(cc, i);
	if (i== cc->index(vh)){
	  make_certificate(f);
	} else if (!has_event(f)){
	  internal::set_facet_label(triangulation(), f, typename Simulator::Event_key());
	}
      }
    }
    v_.new_vertex(vh);
    return vh;
  }

  Cell_handle pop_vertex(Vertex_handle v) {
    v_.delete_vertex(v);
    CGAL_precondition(triangulation_.degree(v)==4);
    std::vector<Cell_handle> ics;
    triangulation_.incident_cells(v, std::back_insert_iterator<std::vector<Cell_handle> >(ics));
    
    for (unsigned int j=0; j< ics.size(); ++j){
      Cell_handle cc= ics[j];
      for (unsigned int i=0; i< 4; ++i){
	Facet f(cc, i);
	if (has_event(f)){
	  simulator()->delete_event(triangulation_.label(f));
	  triangulation_.set_label(f, Event_key());
	}
	for (unsigned int j=0; j<i; ++j){
	  Edge e(cc, i,j);
	  if (has_event(e)){
	    simulator()->delete_event(triangulation_.label(e));
	    triangulation_.set_label(e, Event_key());
	  }
	}
      }
    }
    for (unsigned int j=0; j< ics.size(); ++j){
      Cell_handle cc= ics[j];
      handle_delete_cell(cc);
    }
    v_.delete_cells(ics.begin(), ics.end());
    set_vertex_handle(v->point(), NULL);

    Cell_handle h= triangulation_.tds().remove_from_maximal_dimension_simplex(v);
    handle_new_cell(h);
    v_.new_cells(&h, &h+1);
    for (unsigned int i=0; i< 4; ++i){
      for (unsigned int j=0; j< i; ++j){
	Edge e(h, i, j);
	if (triangulation_.has_degree_3(e)){
	  make_edge_flip(e);
	}
      }
      Facet f(h, i);
      if (!triangulation_.has_degree_3(f)){
	make_certificate(f);
      }
    }
    
    return h;
  }


  void delete_vertex(Point_key){
    CGAL_assertion(0);
  }

  //! dangerous
  /*void suppress_event(const Edge &e) {
    if (triangulation_.label(e) != Simulator::Event_key::null()){
    simulator()->delete_event(triangulation_.label(e) );
    triangulation_.set_label(e, simulator()->null_event());
    }
    }*/


  //! The assertion will catch that the object is in the same sorted order
  /*!
    \todo make this faster.
    This is currently highly suboptimal. Each edge is checked repeatedly if it is degree 3. 
    One posibility would be to track for each edge if it is degree 3 and has a certificate.
  */
  Vertex_handle change_vertex(Point_key k, bool clean=false){
    if (!has_certificates()) return NULL;
    triangulation_.geom_traits().set_time(simulator()->rational_current_time());
    std::vector<Cell_handle> incident_cells;
    triangulation_.incident_cells(vertex_handle(k), back_inserter(incident_cells));
    if (!clean) {
      for (typename std::vector<Cell_handle>::iterator it= incident_cells.begin(); it != incident_cells.end();
	   ++it){
	for (unsigned int i=0; i<4; ++i){
	  Facet f(*it,i);
	  if (!triangulation_.has_degree_3(f) ){
	    Event_key k= triangulation_.label(f);
	    if (has_event(f)) simulator()->delete_event(k);
	    triangulation_.set_label(f, Event_key());
	  }
	  
	  for (unsigned int j=0; j< i; ++j){
	    Edge e(*it, i, j);
	    if (triangulation_.has_degree_3(e) ){
	      Event_key k= triangulation_.label(e);
	      if (has_event(e)) simulator()->delete_event(k);
	      triangulation_.set_label(e, Event_key());
	    }
	  }
	  
	}
      }
      for (typename std::vector<Cell_handle>::iterator it= incident_cells.begin();
	   it != incident_cells.end(); ++it){
	handle_changed_cell(*it);
      }
    } else {
#ifndef NDEBUG
      for (typename std::vector<Cell_handle>::iterator it= incident_cells.begin(); 
	   it != incident_cells.end(); ++it){
	for (unsigned int i=0; i<4; ++i){
	  Facet f(*it,i);
	  CGAL_precondition(!triangulation_.label(f));

	  for (unsigned int j=0; j< i; ++j){
	    Edge e(*it, i, j);
	    CGAL_precondition(!triangulation_.label(e));
	  }
	  
	}
      }
#endif
      for (typename std::vector<Cell_handle>::iterator it= incident_cells.begin();
	   it != incident_cells.end(); ++it){
	handle_new_cell(*it);
      }
      v_.new_cells(incident_cells.begin(), incident_cells.end());
    }

    
    for (typename std::vector<Cell_handle>::iterator it= incident_cells.begin(); it != incident_cells.end();
	 ++it){
      for (unsigned int i=0; i<4; ++i){
	Facet f(*it,i);
	if (!has_event(f)&& !triangulation_.has_degree_3(f)){
	  Event_key k= triangulation_.label(f);
	  make_certificate(f);
	}

	for (unsigned int j=0; j< i; ++j){
	  Edge e(*it, i, j);
	  if ( !has_event(e) && triangulation_.has_degree_3(e)){
	    Event_key k= triangulation_.label(e);
	    make_certificate(e);
	  }
	}
	
      }
    }
    v_.change_vertex(vertex_handle(k));
    return vertex_handle(k);
  }


  typename Triangulation::Vertex_handle new_vertex_regular(Point_key k, Cell_handle h=NULL){
    CGAL_precondition(!has_certificates_);
    typename Triangulation::Vertex_handle vh= triangulation_.insert(k, h);
    v_.new_vertex(vh);
    return vh;
  }


  //! 
  /*!
    Some old certificate edges will be lost, have to find all conflicts first. 

    An important problem case is the edges which were previously degree 3 but are no longer
    after insertion. They can now be degree 4 and will have their certificate deleted properly.
    However, this leaves a face with no certificate (the one or more not destroyed by the new vertex).

    This occurs when a degree 3 edges has a facet added but none destroyed--i.e. a boundary edge
    with degree 3.
  */
  typename Triangulation::Vertex_handle new_vertex(Point_key k, Cell_handle h=NULL){
    /*if (triangulation_.dimension()==3){
      for (All_edges_iterator eit = triangulation_.all_edges_begin();
      eit != triangulation_.all_edges_end(); ++eit){
      if (triangulation_.has_degree_3(*eit)) {
      Event_key k= triangulation_.label(*eit);
      if ( k != Event_key::null()){
      simulator()->delete_event(k);
      triangulation_.set_label(*eit,Event_key::null());
      }
      }
      }
      for (All_facets_iterator eit = triangulation_.all_facets_begin();
      eit != triangulation_.all_facets_end(); ++eit){
      if (!triangulation_.has_degree_3(*eit)) {
      Event_key k= triangulation_.label(*eit);
      if ( k != Event_key::null()){
      simulator()->delete_event(k);
      triangulation_.set_label(*eit,Event_key::null());
      }
      }
      }
      }

      vhs_[k]= triangulation_.insert(k);
      if (triangulation_.dimension()==3){
      initialize_from_scratch();
      }

      return;*/

    Vertex_handle vh;

    if (has_certificates_) {
      triangulation_.geom_traits().set_time(simulator()->rational_current_time());
      if (h==NULL){
	h= triangulation_.locate(k);
      }
      std::vector<Facet> bfacets;
      std::vector<Cell_handle> cells;
      std::vector<Facet> ifacets;
      
      triangulation_.find_conflicts(k, h, back_inserter(bfacets), back_inserter(cells),back_inserter(ifacets));

      /*CGAL_KDS_LOG(LOG_LOTS, "Boundary: ");
      for (unsigned int i=0; i< bfacets.size(); ++i){
	triangulation_.write_labeled_facet(bfacets[i], log_lots() );
      }
      CGAL_KDS_LOG(LOG_LOTS,std::endl);
      CGAL_KDS_LOG(LOG_LOTS,"Internal: ");
      for (unsigned int i=0; i< ifacets.size(); ++i){
	triangulation_.write_labeled_facet(ifacets[i], log_lots() );
      }
      log_lots()  << std::endl;*/

      CGAL_KDS_LOG(LOG_LOTS, "Cells: ");
      for (unsigned int i=0; i< cells.size(); ++i){
	CGAL_KDS_LOG(LOG_LOTS, "[");
	Cell_handle h= cells[i];
	for (unsigned int j=0; j<4; ++j){
	  CGAL_KDS_LOG(LOG_LOTS, h->vertex(j)->point() << " ");
	}
	CGAL_KDS_LOG(LOG_LOTS, "]  ");
      }
      CGAL_KDS_LOG(LOG_LOTS,std::endl);

      //std::sort(ifacets.begin(), ifacets.end());
      //std::sort(bfacets.begin(), bfacets.end());

      /* find boundary edges with degree 3 and propogate their certificate.
	 This may do extra work, so check if the remaining facet is in the interior.
      */
      // propogate it to the facet which will not be destroyed.
      //#if 0
      //std::cout << "ifacets is " << ifacets.size() << " bfacets is " << bfacets.size() << std::endl;
      for (typename std::vector<Facet>::iterator bit = bfacets.begin(); bit != bfacets.end(); ++bit){
	for (unsigned int i=0; i< 3; ++i){
	  //! \todo clean up this mess
	  Edge e= triangulation_.edge(*bit, i);
	  // see if the edges has degree 3
	  if (triangulation_.has_degree_3(e) && has_event(e)){
	    // for simplicity don't try to the push the edge events to facets
	    //simulator()->delete_event(triangulation_.label(e));
	    //triangulation_.set_label(e,  Event_key::null());
	    Facet_circulator fc= triangulation_.incident_facets(e, *bit), fe=fc;
	    ++fc; // to skip the first facet
	    // two facets to look at
	    for (; fc != fe; ++fc){
	      if (!contains(bfacets.begin(), bfacets.end(), *fc)
		  && !contains(ifacets.begin(), ifacets.end(), *fc)){
		// we have a facet to repair which is not obsolete, and not on the boundary
		// must check if another edge is degree 3
		bool d3=false;
		for (int k=0; k<3; ++k){
		  Edge e2= triangulation_.edge(*fc, k);
		  // have to make sure that the other edge has not gone away with this
		  if (triangulation_.equal(e, e2)) continue;
		  else if ( triangulation_.has_degree_3(e2)
			   && triangulation_.label(e2) ){
		    d3= true;
		    break;
		  }
		}
		if (!d3){
		  CGAL_KDS_LOG(LOG_LOTS,"Handled lost edge: was ");
		  //triangulation_.write_labeled_edge(e, log_lots());
		  CGAL_KDS_LOG(LOG_LOTS, " moves to facet ");
		  //triangulation_.write_facet(*fc, log_lots());
		  CGAL_KDS_LOG(LOG_LOTS,std::endl);
		  
		  typename Simulator::Event_key key= triangulation_.label(e);
		  Event_key k= change_to_facet_flip(*fc, key);

		  CGAL_KDS_LOG(LOG_LOTS, "Event is " << k << std::endl);
		  
		  //triangulation_.set_label(e, P::Event_key::null());
		  //typename P::Event_key k=simulator()->new_event();
		  triangulation_.set_label(*fc, k);
		  
		  triangulation_.set_label(e, typename Simulator::Event_key());
		} else {
		  typename Simulator::Event_key k= triangulation_.label(e);
		  simulator()->delete_event(k);
		  triangulation_.set_label(e, typename Simulator::Event_key());
		}
		break; // only can have one to deal with
	      }
	    }
	  }
	}
      }
      //#endif
      
 
      /* have to delete all facet and edge certificates
	 Could use the other lists, but the easiest for the edges is to just use the cells
      */
      for (typename std::vector<Cell_handle>::iterator cit= cells.begin(); cit != cells.end(); ++cit){
	for (unsigned int i=0; i< 4; ++i){
	  if (has_event(Facet(*cit, i))){
	    typename Simulator::Event_key k= triangulation_.label(Facet(*cit,i));
	    simulator()->delete_event(k);
	    triangulation_.set_label(Facet(*cit, i), typename Simulator::Event_key());
	  }
	  for (unsigned int j=0; j<i; ++j){
	    Edge e(*cit, i,j);
	    if (has_event(e)){
	      typename Simulator::Event_key k= triangulation_.label(e);
	      simulator()->delete_event(k);
	      triangulation_.set_label(e, typename Simulator::Event_key());
	    }
	  }
	}
      }

      for (typename std::vector<Cell_handle>::iterator cit= cells.begin(); cit != cells.end(); ++cit){
	handle_delete_cell(*cit);
      }
      v_.delete_cells(cells.begin(), cells.end());
      //! \todo replace by insert_in_hole
      CGAL_KDS_LOG(LOG_LOTS, "Inserting.\n");
      vh=triangulation_.insert(k, h);
      set_vertex_handle(k,vh);
      change_vertex(k, true);
    } else {
      // if low dimensional
      vh=triangulation_.insert(k, h);
    }
    
    CGAL_postcondition(audit_structure());
    v_.new_vertex(vh);
    return vh;
  }



  bool has_certificates() const {
    return has_certificates_;
  }

  void set_has_certificates(bool b) {
    if (has_certificates_ && !b){
      destroy_all_certificates();
     } else if (!has_certificates_ && b){
       if (triangulation().dimension() == 3){
	create_all_certificates();
      }
    }
  }













  Facet flip(const Edge &edge){
    v_.pre_edge_flip(edge);
    if (true) {
      CGAL_KDS_LOG(LOG_LOTS,"\n\nFlipping edge ");
      //triangulation_.write_labeled_edge(edge, log_lots());
      CGAL_KDS_LOG(LOG_LOTS,std::endl);
    }
    CGAL_assertion(triangulation_.tds().is_edge(edge.first, edge.second, edge.third) || print());
            
    Vertex_handle poles[2];
    poles[0]= edge.first->vertex(edge.second);
    poles[1]= edge.first->vertex(edge.third);
    int polesi[2];
    polesi[0] = edge.first->index(poles[0]);
    polesi[1]= edge.first->index(poles[1]);

    typename Simulator::Event_key failed_key = triangulation_.label(edge);
    CGAL_precondition(failed_key);
    typename Simulator::Root_stack ore= extract_root_stack(failed_key);
    
    Facet neighboring_facet= triangulation_.opposite(Facet(edge.first, polesi[0]));
    
    triangulation_.set_label(edge, typename Simulator::Event_key());

    // handle the cross edges to make sure that they are no longer edge flips
    {
      Cell_circulator cc= triangulation_.incident_cells(edge), ce= cc;
      do {
	Cell_handle cell=cc;
	Edge eic= triangulation_.edge_in_cell(edge, cell);
	Edge xe= triangulation_.cross(eic);
	if (triangulation_.has_degree_3(xe)){
	  make_not_edge_flip(xe, cc);
	}
	// Just try them all
	for (unsigned int i=0; i<4; ++i){
	  Facet f(cc, i);
	  if (has_event(f)){
	    simulator()->delete_event(triangulation_.label(f));
	    triangulation_.set_label(f, typename Simulator::Event_key());
	  }

	  // hack for recently redundant vertices
	  for (unsigned int j=0; j<i; ++j){
	    Edge e(cell, i,j);
	    if (has_event(e)){
	      simulator()->delete_event(triangulation_.label(e));
	      triangulation_.set_label(e, typename Simulator::Event_key());
	    }
	  }
	}
      } while (++cc != ce);
    }
#ifndef NDEBUG
    {
      Cell_circulator cc= triangulation_.incident_cells(edge), ce= cc;
      do {
	Cell_handle cell=cc;
	for (unsigned int i=0; i<4; ++i){
	  Facet f(cell, i);
	  CGAL_assertion(!has_event(f));
	  for (unsigned int j=0; j<i; ++j){
	    Edge e(cell, i,j);
	    if (has_event(e)){
	      triangulation_.write_labeled_edge(e, std::cerr);
	    }
	    CGAL_assertion(!has_event(e));
	  }
	}
      } while (++cc != ce);
    }
#endif
    {
      Cell_circulator cc= triangulation_.incident_cells(edge), ce= cc;
      do {
	handle_delete_cell(cc);
	v_.delete_cells(&cc, &cc+1);
      } while (++cc != ce);
      
    }

    triangulation_.tds().flip_flippable(edge);
    CGAL_expensive_assertion(labeling_is_valid());
    /*if (verbose){
      write_labeled_state();
      }*/

    CGAL_assertion(triangulation_.tds().is_facet(neighboring_facet.first, neighboring_facet.second));
    Facet internal_facet= triangulation_.opposite(neighboring_facet);

   
      
    Cell_handle cells[2];
    cells[1]= internal_facet.first;
    polesi[1]= cells[1]->index(poles[1]);
    cells[0]= cells[1]->neighbor(polesi[1]);
    polesi[0]= cells[0]->index(poles[0]);
      
    Facet middle_facet(cells[0], polesi[0]);
    triangulation_.clear_cell_labels(middle_facet.first);
    triangulation_.clear_cell_labels(triangulation_.opposite(middle_facet).first);

    if (true){
      CGAL_KDS_LOG(LOG_LOTS, "Created facet ");
      //triangulation_.write_labeled_facet(Facet(cells[0], polesi[0]), log_lots());
      CGAL_KDS_LOG(LOG_LOTS, std::endl);
    }
    
    for (int c=0; c<2; ++c){
      handle_new_cell(cells[c]);
    }
    v_.new_cells(&cells[0], &cells[0]+2);

    bool id3[2][3];
    // do certifcates for d3 edges
    for (int c=0; c<2; ++c){
      for (int i=0; i<3; ++i){
	Edge e= triangulation_.edge_around_vertex(cells[c], poles[c], i);
	if (triangulation_.has_degree_3(e)){
	  make_edge_flip(e);
	  id3[c][i]=true;
	} else {
	  id3[c][i]=false;
	}
      }
    }
    
    // build certificates for faces
    //! \todo make this more efficient
    /*! 
      Edges are checked too many times to see if they are degree 3
    */
    for (int c=0; c<2; ++c){
      for (int i=0; i<3; ++i){
	Facet f= triangulation_.facet_around_vertex(cells[c], poles[c], i);
	CGAL_assertion(f.second != cells[c]->index(poles[c]));
	bool hd3= (id3[c][(i+1)%3] || id3[c][(i+2)%3]);
	CGAL_assertion(hd3== triangulation_.has_degree_3(f));
	if (!hd3) make_certificate(f);
      }
    }


    if (!ore.empty()){
      typename Simulator::Time t= ore.top();
      ore.pop();
      typename Simulator::Event_key k= simulator()->new_event(t, Facet_flip_event(ore, middle_facet, tr_.wrapper_pointer()));

      triangulation_.set_label(middle_facet, k);
    } else {
      triangulation_.set_label(middle_facet, simulator()->null_event());
    }
    //generate_geometry();

    CGAL_postcondition(audit_structure());
    v_.post_edge_flip(middle_facet);
    return middle_facet;
  }
















  Edge flip(const Facet &flip_facet){
    v_.pre_facet_flip(flip_facet);
    Vertex_handle poles[2];
    Facet other_flip_facet= triangulation_.opposite(flip_facet);
    poles[0]= flip_facet.first->vertex(flip_facet.second);
    poles[1]= other_flip_facet.first->vertex(other_flip_facet.second);
      
    if (true) {
      CGAL_KDS_LOG(LOG_LOTS,"\n\nFlipping facet ");
      //triangulation_.write_labeled_facet(flip_facet, log_lots() );
      CGAL_KDS_LOG(LOG_LOTS," with poles " << poles[0]->point() << ", " << poles[1]->point());
      CGAL_KDS_LOG(LOG_LOTS,std::endl);
    }
    CGAL_assertion(triangulation_.tds().is_facet(flip_facet.first, flip_facet.second) || print());

    Event_key failed_key= triangulation_.label(flip_facet);

    typename Simulator::Root_stack ore= extract_root_stack(failed_key);
    CGAL_KDS_LOG(LOG_LOTS, "Next root of solver is " << ore.top() << std::endl);

    //typename P::Event_key index= triangulation_.label(flip_facet);

    Facet inside_facet(flip_facet.first, (flip_facet.second+1)%4);
    Facet neighbor_facet= triangulation_.opposite(inside_facet);
    Cell_handle cells[2]= {flip_facet.first, other_flip_facet.first};
    CGAL_assertion(neighbor_facet.first != cells[0] && neighbor_facet.first != cells[1]);
     
    // go around and change the handler for each edge if it is degree 3. Don't have to look off cell
    for (int c=0; c<2; ++c){
      for (int i=0; i<3; ++i){
	Edge e= triangulation_.edge_around_vertex(cells[c], poles[c], i);
	if (triangulation_.has_degree_3(e)){
	  make_not_edge_flip(e, cells[c]);
	}

	Facet f = triangulation_.facet_around_vertex(cells[c], poles[c],i);
	if (has_event(f)){
	  Event_key k= triangulation_.label(f);
	  CGAL_KDS_LOG(LOG_LOTS,"Cleaning facet ");
	  //triangulation_.write_labeled_facet(f, log_lots() );
	  CGAL_KDS_LOG(LOG_LOTS,std::endl);
	  triangulation_.set_label(f, typename Simulator::Event_key());
	  simulator()->delete_event(k);
	}
      }
    }
    for (int c=0; c<2; ++c){
      handle_delete_cell(cells[c]);
    }
    v_.delete_cells(&cells[0], &cells[0]+2);
#ifndef NDEBUG
    CGAL_KDS_LOG(LOG_LOTS, *this);
#endif
    triangulation_.tds().flip_flippable(flip_facet);

    CGAL_assertion(triangulation_.tds().is_facet(neighbor_facet.first, neighbor_facet.second));

    Facet a_facet= triangulation_.opposite(neighbor_facet);
    Cell_handle a_cell= a_facet.first;

#ifndef NDEBUG
    if (!a_cell->has_vertex(poles[0]) || !a_cell->has_vertex(poles[1])){
      CGAL_KDS_LOG(LOG_LOTS, "Cell: ");
      for (int i=0; i<4; ++i){
	CGAL_KDS_LOG(LOG_LOTS, a_cell->vertex(i)->point() << " ");
      }
      //triangulation_.write_labeled_facet(a_facet, log_lots());
      CGAL_KDS_LOG(LOG_LOTS,std::endl);
      //triangulation_.write_labeled_facet(neighbor_facet, log_lots());
      CGAL_KDS_LOG(LOG_LOTS, std::endl);
      CGAL_assertion_code(print());
      CGAL_assertion(0);
    }
#endif

    Edge central_edge(a_cell, a_cell->index(poles[0]), a_cell->index(poles[1]));

    Cell_circulator cc= triangulation_.incident_cells(central_edge), ce=cc;
    do {
      triangulation_.clear_cell_labels(cc);
    } while (++cc != ce);
#ifndef NDEBUG
    CGAL_KDS_LOG(LOG_LOTS, *this);
#endif
    CGAL_expensive_assertion(labeling_is_valid());

    if (true) {
      CGAL_KDS_LOG(LOG_LOTS, "Created edge ");
      //triangulation_.write_labeled_edge(central_edge, log_lots() );
      CGAL_KDS_LOG(LOG_LOTS,std::endl);
    }
    {
      Cell_circulator cc= triangulation_.incident_cells(central_edge), ce=cc;
      do {
	handle_new_cell(cc);
	v_.new_cells(&cc, &cc+1);
      } while (++cc != ce);
    }
    // fix the cross edges, must make them non facet flip if necessary
    {
      Cell_circulator cc= triangulation_.incident_cells(central_edge), ce=cc;
      do {
	Edge eic= triangulation_.edge_in_cell(central_edge, cc);
	Edge xe= triangulation_.cross(eic);
	bool d3= triangulation_.has_degree_3(xe);
	if (d3){
	  make_no_events(xe);
	  make_certificate(xe);
	} else {
	  // check adjacent facets
	  make_certificate(Facet(cc, eic.second));
	  make_certificate(Facet(cc, eic.third));
	}
      } while (++cc != ce);
    }

    //typename Simulator::Time t= s.next_time_negative();
    if (!ore.empty()){
      typename Simulator::Time t= ore.top();
      ore.pop();
      typename Simulator::Event_key k= simulator()->new_event(t, Edge_flip_event(ore, central_edge, tr_.wrapper_pointer()));
      triangulation_.set_label(central_edge, k);
    } else {
      triangulation_.set_label(central_edge, simulator()->null_event());
    }

    //generate_geometry();

    CGAL_postcondition(audit_structure());
    v_.post_facet_flip(central_edge);
    return central_edge;
  }










  void audit() const {
    CGAL_KDS_LOG(LOG_SOME, "Verifying at time " << simulator()->audit_time() << ".\n");
    triangulation_.geom_traits().set_time(simulator()->rational_current_time());
    //CGAL_precondition(triangulation_.geom_traits().time()== simulator()->audit_time());
    //triangulation_.geom_traits().set_time(simulator().rational_current_time());
    audit_structure();
    triangulation_.is_valid(true);
  }









  
  const Triangulation& triangulation() const {
    return triangulation_;
  }

  Triangulation& triangulation() {
    return triangulation_;
  }







 

  const Moving_object_table* moving_object_table() const {
    return tr_.moving_point_table_pointer();
  }

  Simulator* simulator() {
    return tr_.simulator_pointer();
  }

  const Simulator* simulator() const {
    return tr_.simulator_pointer();
  }
  
  const Point& point(Point_key k) const { 
    return moving_object_table()->at(k);
  }

  Vertex_handle vertex_handle(Point_key k) const {
    if (k.index()>=0 && static_cast<unsigned int>(k.index()) < vhs_.size()){
      return vhs_[k.index()];
    } else {
      return NULL;
    }
  }
  void set_vertex_handle(Point_key k, Vertex_handle vh) {
    CGAL_precondition(k.index() >=0);
    unsigned int bin= k.index();
    while (vhs_.size() <= bin){
      vhs_.push_back(NULL);
    }
    /*if (vhs_.size() <=bin){
      vhs_.resize(bin+1);
      }*/
    vhs_[k.index()]=vh;
  }

  typename TraitsT::Side_of_oriented_sphere_3 power_test_object() const {
    return soc_;
  };
  typename TraitsT::Orientation_3 orientation_object() const {
    return o3_;
  }
  typename Simulator::Root_stack extract_root_stack(Event_key k) const {
    //typename Simulator::Event_handle<Event_base> eh(simulator()->event(k, Event_base()));
    //typename Simulator::Root_stack s= eh.pointer()->root_stack();
    
    return simulator()->event(k, Event_base()).root_stack();
  }
  /*
  typename Simulator::Time extract_time(Event_key k) const {
    return simulator()->event<Facet_flip_event::Base>(k)->time();
    }*/
  
private:
  typename Simulator::Root_stack root_stack(const Edge &e) const {
    return root_stack(*triangulation_.incident_facets(e));
  }


  void destroy_all_certificates() {
    //vhs_.clear();
    CGAL_precondition(has_certificates_);
    for (All_edges_iterator eit = triangulation_.all_edges_begin();
	 eit != triangulation_.all_edges_end(); ++eit){
      Event_key k= triangulation_.label(*eit);
      if ( k ){
	simulator()->delete_event(k);
	triangulation_.set_label(*eit,Event_key());
      }
    }
    for (All_facets_iterator eit = triangulation_.all_facets_begin();
	 eit != triangulation_.all_facets_end(); ++eit){
      Event_key k= triangulation_.label(*eit);
      if ( k ){
	simulator()->delete_event(k);
	triangulation_.set_label(*eit,Event_key());
      }
      //}
    }
    for (All_cells_iterator eit = triangulation_.all_cells_begin();
	 eit != triangulation_.all_cells_end(); ++eit){
      handle_delete_cell(eit);
    }
    v_.delete_cells(triangulation_.all_cells_begin(), triangulation_.all_cells_end());
    has_certificates_=false;
  }









  void create_all_certificates() {
    CGAL_precondition(!has_certificates_);
    for (All_cells_iterator eit = triangulation_.all_cells_begin();
	 eit != triangulation_.all_cells_end(); ++eit){
      handle_new_cell(eit);
    }
    v_.new_cells(triangulation_.all_cells_begin(), triangulation_.all_cells_end());
    for (All_edges_iterator eit = triangulation_.all_edges_begin(); 
	 eit != triangulation_.all_edges_end(); ++eit){
      if (triangulation_.has_degree_3(*eit)) {
	CGAL_KDS_LOG(LOG_LOTS,"initializing [");
	//triangulation_.write_labeled_edge(*eit, log_lots() );
	CGAL_KDS_LOG(LOG_LOTS, "]\n");
	make_certificate(*eit);
      }
    }
    for (All_facets_iterator eit = triangulation_.all_facets_begin(); 
	 eit != triangulation_.all_facets_end(); ++eit){
      if (!triangulation_.has_degree_3(*eit)) {
	CGAL_KDS_LOG(LOG_LOTS, "initializing [");
	//triangulation_.write_labeled_facet(*eit, log_lots() );
	CGAL_KDS_LOG(LOG_LOTS, "]\n");
	make_certificate(*eit);
      }
    }
    //CGAL_precondition(vhs_.empty());
    for (Finite_vertices_iterator vit = triangulation_.finite_vertices_begin(); vit != triangulation_.finite_vertices_end(); ++vit){
      set_vertex_handle(vit->point(), vit);
    }
    has_certificates_=true;
  }





  void make_no_events(const Edge &e) {
    CGAL_precondition(triangulation_.has_degree_3(e));
    Facet_circulator fc= triangulation_.incident_facets(e);
    Facet_circulator fe= fc;
    do {
      if (has_event(*fc)){
	Event_key k= triangulation_.label(*fc);
	simulator()->delete_event(k);
	triangulation_.set_label(*fc, Event_key());
      }
    } while(++fc != fe);
  }






  void make_certificate( const Edge &e){
    CGAL_KDS_LOG(LOG_LOTS, "makeing certificate for edge ");
    CGAL_KDS_LOG_WRITE(LOG_LOTS, triangulation_.write_edge(e, LOG_STREAM));
    CGAL_KDS_LOG(LOG_LOTS, std::endl);
    CGAL_precondition(triangulation_.has_degree_3(e));
    CGAL_precondition_code(Facet_circulator fc= triangulation_.incident_facets(e));
    CGAL_precondition_code(Facet_circulator fe= fc);
    CGAL_precondition_code(do {);
			   CGAL_precondition(!has_event(*fc));
			   CGAL_precondition_code(++fc);
			   CGAL_precondition_code(}while(fc != fe));
    //! one endpoint is the center of an isolated tet, so this edge can't flip anyway.
    /*if (triangulation_.degree(triangulation_.vertex(e,1))==4 
      || triangulation_.degree(triangulation_.vertex(e,0))==4 ){
      triangulation_.set_label(e, simulator()->null_event());
      return;
      }*/
    typename Simulator::Root_stack s= root_stack(e);
    //typename Simulator::Time t= s.next_time_negative();
    if (!s.empty()){
      CGAL_KDS_LOG(LOG_LOTS,"Failure time is " << s.top() << std::endl);
      typename Simulator::Time t= s.top();
      s.pop();
      CGAL_KDS_LOG(LOG_LOTS, "Next root of this cert is " << s.top() << std::endl);
      typename Simulator::Event_key k= simulator()->new_event(t, Edge_flip_event(s, e, tr_.wrapper_pointer()));
      triangulation_.set_label(e, k);
    } else {
      triangulation_.set_label(e, simulator()->null_event());
    }
  }

  bool has_event(const Edge &e) const {
    return triangulation_.label(e);
  }

  bool has_event(const Facet &e) const {
    return triangulation_.label(e);
  }

  void make_certificate( const Facet &e){
    CGAL_KDS_LOG(LOG_LOTS, "makeing certificate for facet ");
    //triangulation_.write_facet(e, log_lots());
    CGAL_KDS_LOG(LOG_LOTS,  std::endl);


    CGAL_precondition(!triangulation_.has_degree_3(e));
    for (int i=0; i<3; ++i){
      CGAL_precondition(!triangulation_.label(triangulation_.edge(e, i)));
    }
    typename Simulator::Root_stack s= root_stack(e);
    //typename Simulator::Time t= s.next_time_negative();
    if (!s.empty()){
      typename Simulator::Time t= s.top();
      s.pop();
      CGAL_KDS_LOG(LOG_LOTS, "Next root of this cert is " << s.top() << std::endl);
      typename Simulator::Event_key k= simulator()->new_event(t, Facet_flip_event(s, e, tr_.wrapper_pointer()));
      triangulation_.set_label(e, k);
    } else {
      triangulation_.set_label(e, simulator()->null_event());
    }
  }

  
  template <class Oit>
  void point_keys(const Facet &f, Oit out) const {
    int hinf=-1;
    for (unsigned int i=0; i<4; ++i){
      Point_key k= f.first->vertex(i)->point();
      if (!k){
	hinf=i;
	break;
      }
    }
    if (hinf==-1){
      Point_key k= triangulation_.mirror_vertex(f.first, f.second)->point();
      if ( !k ){
	hinf=4;
      }
    }
    if (hinf != -1) {
      CGAL_KDS_LOG(LOG_LOTS, "hinf is " << hinf << std::endl);
      if (hinf ==4){
	for (unsigned int i=0; i<4; ++i){
	  Point_key k= f.first->vertex(i)->point();
	  *out= k;
	  ++out;
	}
	return;
      } else {
	//Facet ff(f.first, hinf);
	if (hinf%2!=0){
	  Point_key k= triangulation_.mirror_vertex(f.first, f.second)->point();
	  *out= k;
	  ++out;
	}
	for (int i=0; i<4; ++i){
	  // CGAL infinite cells seem to be misoriented
	  if (i==hinf) continue;
	  *out= f.first->vertex(i)->point();
	  ++out;
	}
	if (hinf%2==0){
	  Point_key k= triangulation_.mirror_vertex(f.first, f.second)->point();
	  *out= k;
	  ++out;
	}
      }
    } else {
      for (unsigned int i=0; i<4; ++i){
	Point_key k= f.first->vertex(i)->point();
	*out= k;
	++out;
      }
      Point_key k= triangulation_.mirror_vertex(f.first, f.second)->point();
      *out= k;
      ++out;
    }
  }

 

  void make_edge_flip(Edge &edge){
    CGAL_KDS_LOG(LOG_LOTS, "Making edge flip ");
    //triangulation_.write_edge(edge, log_lots() );
    CGAL_KDS_LOG(LOG_LOTS,std::endl);
    CGAL_assertion(triangulation_.has_degree_3(edge));
    typename Simulator::Event_key k= typename Simulator::Event_key();
    Facet_circulator fc= triangulation_.incident_facets(edge), fe=fc;
    do {
      if (has_event(*fc)){
	CGAL_assertion( !k );
	k=triangulation_.label(*fc);
	triangulation_.set_label(*fc, typename Simulator::Event_key());
#ifdef NDEBUG
	break;
#endif
      }
    } while (++fc != fe);

    //CGAL_assertion(k!= Event_key::null());
    /*if (k== Simulator::Event_key::null()){
      log_lots()  << "Non failing certificate on edge ";
      triangulation_.write_labeled_edge(edge, log_lots() );
      log_lots()  << std::endl;
      }*/
    if ( k ){
      /*if (triangulation_.degree(triangulation_.vertex(edge,0))==4 
	|| triangulation_.degree(triangulation_.vertex(edge,1))==4){
	triangulation_.set_label(edge, simulator()->null_event());
	} else {*/
      Event_key kk=change_to_edge_flip(edge, k);
      triangulation_.set_label(edge, kk);
      //}
    } else {
      CGAL_KDS_LOG(LOG_LOTS, "Making up edge event.\n");
      make_certificate(edge);
    }
  }

  void make_not_edge_flip(Edge &edge, Cell_handle cell){
    if (true){
      CGAL_KDS_LOG(LOG_LOTS, "Making edge ");
      //triangulation_.write_labeled_edge(edge, log_lots() );
      CGAL_KDS_LOG(LOG_LOTS, " not an edge flip.\n");
    }
    CGAL_assertion(triangulation_.has_degree_3(edge) || print());
    typename Simulator::Event_key index = triangulation_.label(edge);
    CGAL_precondition( index );

    triangulation_.set_label(edge, typename Simulator::Event_key());
    
    Facet_circulator fc=triangulation_.incident_facets(edge), ef= fc;
    do {
      Facet cf= *fc;
      if (cf.first == cell || triangulation_.opposite(cf).first == cell){
	CGAL_assertion( !triangulation_.label(cf) );
      } else {
	bool hd3=false;
	for (int i=0; i< 3; ++i){
	  Edge ec= triangulation_.edge(cf, i);
	  if (triangulation_.equal(ec, edge)) continue;
	  else {
	    if (triangulation_.has_degree_3(ec)){
	      hd3=true;
	      break;
	    }
	  }
	}
	
	if (!hd3){
	  Event_key k= change_to_facet_flip(cf, index);
	  triangulation_.set_label(cf, k);
	  CGAL_KDS_LOG(LOG_LOTS,"put cert on ");
	  CGAL_KDS_LOG_WRITE(LOG_LOTS, triangulation_.write_facet(cf,LOG_STREAM  ));
	  CGAL_KDS_LOG(LOG_LOTS, std::endl);
	} else {
	  simulator()->delete_event(index);
	  if (true){
	    CGAL_KDS_LOG(LOG_LOTS, "Facet ");
	    CGAL_KDS_LOG_WRITE(LOG_LOTS, triangulation_.write_labeled_facet(cf, LOG_STREAM) );
	    CGAL_KDS_LOG(LOG_LOTS, " has another degree three edge, deleting certificate.\n");
	  }
	}
      }
    } while (ef != ++fc);
  }

  template <class C>
  std::back_insert_iterator<C> back_inserter(C &c) const {
    return std::back_insert_iterator<C>(c);
  }
  
  bool contains(typename std::vector<Facet>::iterator beg, 
		typename std::vector<Facet>::iterator end, Facet f){
    Facet of = triangulation_.opposite(f);
    for (; beg != end; ++beg){
      if (f.first == beg->first){
	if (f.second == beg->second) return true;
      } else if (of.first == beg->first){
	if (of.second == beg->second) return true;
      }
    }
    return false;
  }

  bool labeling_is_valid() const {
    for (All_facets_iterator eit= triangulation_.all_facets_begin(); eit != triangulation_.all_facets_end();
	 ++eit){
      Facet f= *eit;
      Facet of= triangulation_.opposite(f);
      if (triangulation_.label(f) != triangulation_.label(of)){
	triangulation_.write_labeled_facet(f, std::cerr);
	std::cerr << " does not match ("<< triangulation_.label(f) << ", " 
		  << triangulation_.label(of) << ")\n";
      }
    }
    for (All_edges_iterator eit= triangulation_.all_edges_begin(); eit != triangulation_.all_edges_end();
	 ++eit){
      Edge e=*eit;
      triangulation_.label(e);
    }
    return true;
  }

  bool audit_structure() const {
    if (triangulation_.dimension() != 3) return true;
    //print();
    
    for (All_edges_iterator eit = triangulation_.all_edges_begin(); 
	 eit != triangulation_.all_edges_end(); ++eit){
      if (!triangulation_.has_degree_3(*eit)){
	if (has_event(*eit)){
	  triangulation_.write_labeled_edge(*eit, std::cerr);
	}
	CGAL_assertion(!has_event(*eit));
      } else {
	if (has_certificates_ && !has_event(*eit)){
	  triangulation_.write_labeled_edge(*eit, std::cerr);
	}
	CGAL_assertion(!has_certificates_ || has_event(*eit));
	//CGAL_assertion(triangulation_.label(*eit)!= Simulator::Event_key::null());
      }
    }
    for (All_facets_iterator eit = triangulation_.all_facets_begin(); 
	 eit != triangulation_.all_facets_end(); ++eit){
      if (triangulation_.has_degree_3(*eit)){
	if (has_event(*eit)){
	  triangulation_.write_labeled_facet(*eit, std::cerr);
	}
	CGAL_assertion(!has_event(*eit));
      } else {
	if (has_certificates_ && !has_event(*eit)){
	  triangulation_.write_labeled_facet(*eit, std::cerr);
	}
	CGAL_assertion(!has_certificates_|| has_event(*eit));
      }
    }
    return true;
  }


  Event_key change_to_edge_flip(const Edge &e, Event_key k){
    if (k== simulator()->null_event()) return k;
    typename Simulator::Root_stack s= extract_root_stack(k);
    return simulator()->set_event(k, Edge_flip_event(s, e, tr_.wrapper_pointer()));
  }

  Event_key change_to_facet_flip(const Facet &f, Event_key k){
    if (k== simulator()->null_event()) return k;
    typename Simulator::Root_stack s= extract_root_stack(k);
    return simulator()->set_event(k, Facet_flip_event(s, f, tr_.wrapper_pointer()));
  }

  /*Moving_object_table* moving_object_table() {
    return mpt_.pointer();
    }*/

#ifdef DELAUNAY_CACHING
  void handle_delete_cell(Cell_handle h){
    h->info().clear();
  }

  void handle_new_cell(Cell_handle){
    //h->info().update(h, moving_object_table());
  }
  void handle_changed_cell(Cell_handle h){
    h->info().clear();
    handle_new_cell(h);
  }

  struct Full_certificate {
    Full_certificate(Cell_handle h, Point pt): h_(h), pt_(pt){}
    template <class CC>
    typename CC::result_type operator()(const CC &cc) const {
   
      CGAL_KDS_LOG_MAPLE( std::endl << std::endl);
      CGAL_KDS_LOG_MAPLE("(" << cc(CGAL::KDS::internal::point(pt_).x())<< ")*(" << h_->info().x());
      CGAL_KDS_LOG_MAPLE(")+ (" << cc(CGAL::KDS::internal::point(pt_).y()) << ")*(" << h_->info().y());
      CGAL_KDS_LOG_MAPLE(")+ (" << cc(CGAL::KDS::internal::point(pt_).z())<< ")*("<< h_->info().z());
      CGAL_KDS_LOG_MAPLE(")+(" <<CGAL::KDS::internal::lift(pt_, cc)<<")*("<<h_->info().l());
      CGAL_KDS_LOG_MAPLE(")+(" << h_->info().o() << ");\n");
      return
	cc(CGAL::KDS::internal::point(pt_).x())*h_->info().x()
	+cc(CGAL::KDS::internal::point(pt_).y())*h_->info().y()
	+cc(CGAL::KDS::internal::point(pt_).z())*h_->info().z()
	+CGAL::KDS::internal::lift(pt_, cc)*h_->info().l()
	+h_->info().o();
    }
    Cell_handle h_;
    Point pt_;
  };

  struct Hull_certificate {
    Hull_certificate(Cell_handle h, Point pt): h_(h), pt_(pt){
    }
    template <class CC>
    typename CC::result_type operator()(const CC &cc) const {
      CGAL_KDS_LOG_MAPLE(std::endl << std::endl);
      CGAL_KDS_LOG_MAPLE("("<< cc(CGAL::KDS::internal::point(pt_).x()) << ")*(" << h_->info().x());
      CGAL_KDS_LOG_MAPLE(")+(" << cc(CGAL::KDS::internal::point(pt_).y())<<")*("<<h_->info().y());
      CGAL_KDS_LOG_MAPLE(")+("<<cc(CGAL::KDS::internal::point(pt_).z())<<")*("<<h_->info().z());
      CGAL_KDS_LOG_MAPLE(")+("<<h_->info().o()<<");\n");
      return 
	cc(CGAL::KDS::internal::point(pt_).x())*h_->info().x()
	+cc(CGAL::KDS::internal::point(pt_).y())*h_->info().y()
	+cc(CGAL::KDS::internal::point(pt_).z())*h_->info().z()
	+h_->info().o();
    }
    Cell_handle h_;
    Point pt_;
  };
  
  typename Simulator::Root_stack root_stack(const Facet &f) const {
    typename Simulator::Root_stack s;
#ifndef NDEBUG
    root_stack_nc(f);
    
#endif
    if ( !f.first->vertex(f.second)->point() ){
      // this is inf
      Point_key mv= triangulation_.mirror_vertex(f.first, f.second)->point();
      if (!f.first->info().set()){
	f.first->info().update(f.first, moving_object_table());
      }
      s= simulator()->root_stack_object(Hull_certificate(f.first,
						     moving_object_table()->at(mv)));
    } else {
      bool hinf=false;
      Cell_handle nc= f.first->neighbor(f.second);
      for (unsigned int i=0; i<4; ++i){
	if ( !nc->vertex(i)->point() ){
	  hinf=true;
	  break;
	}
      }

      
      if (hinf){
	Point_key mv= f.first->vertex(f.second)->point();
	Cell_handle h= f.first->neighbor(f.second);
	if (!h->info().set()){
	  h->info().update(h, moving_object_table());
	}

	s= simulator()->root_stack_object(Hull_certificate(h,
						       moving_object_table()->at(mv)));
      } else {
	Cell_handle h= f.first;
	if (h->info().set()){
	  Point_key mv= triangulation_.mirror_vertex(f.first, f.second)->point();

	  s= simulator()->root_stack_object(Full_certificate(h,
							 moving_object_table()->at(mv)));
	} else {
	  Cell_handle h= f.first->neighbor(f.second);
	  if (!h->info().set()){
	    h->info().update(h, moving_object_table());
	  }
	  Vertex_handle v= f.first->vertex(f.second);




	  s= simulator()->root_stack_object(Full_certificate(h,
							 moving_object_table()->at(v->point())));
	}
      }
    }
    
    //typename Simulator::Root_stack snc= root_stack_nc(f);
    
    return s;
    
  }
#else
  void handle_delete_cell(Cell_handle){
  }
  
  void handle_new_cell(Cell_handle){
  }
  void handle_changed_cell(Cell_handle){

  }

  

  typename Simulator::Root_stack root_stack(const Facet &f) const {
    return root_stack_nc(f);
  }
#endif

  typename Simulator::Root_stack root_stack_nc(const Facet &f) const {
    std::vector<Point_key> ids;
    point_keys(f, back_inserter(ids));
#ifndef NDEBUG
    std::vector<Point_key> mids;
    Facet of= triangulation_.opposite(f);
    point_keys(of, back_inserter(mids));
#endif
#ifndef NDEBUG
    CGAL_KDS_LOG(LOG_LOTS, "Creating root_stack for points ");
    for (typename std::vector<Point_key>::const_iterator cit= ids.begin(); cit != ids.end(); ++cit){
      CGAL_KDS_LOG(LOG_LOTS, *cit);
    }
    CGAL_KDS_LOG(LOG_LOTS, std::endl);
#endif
    bool is_const=true;
    for (unsigned int i=0; i<ids.size(); ++i){
      if (!point(ids[i]).is_constant()){
	is_const=false;
	break;
      }
    }
    if (is_const) {
      typename Kinetic_kernel::Certificate_function cf(1.0);
      return simulator()->root_stack_object(cf);
    } else if (ids.size()==4){
      /*if (point(ids[0]).is_constant()){
      // hack
      std::swap(ids[0], ids[3]);
      std::swap(ids[1], ids[2]);
      }*/
      return simulator()->root_stack_object(o3_(point(ids[0]),
						     point(ids[1]),
						     point(ids[2]),
						     point(ids[3])) );
    } else {
      CGAL_assertion(ids.size()==5);
      /*if (point(ids[0]).is_constant()){
      // hack for linear case
      std::swap(ids[0], ids[4]);
      std::swap(ids[1], ids[2]);
      }*/
      return simulator()->root_stack_object(soc_(point(ids[0]),
						      point(ids[1]),
						      point(ids[2]),
						      point(ids[3]),
						      point(ids[4])));
    }
    
  }

  TraitsT tr_;
  Triangulation triangulation_;
  std::vector<Vertex_handle> vhs_;
  typename TraitsT::Side_of_oriented_sphere_3 soc_;
  typename TraitsT::Orientation_3 o3_;
  bool has_certificates_;
  Visitor v_;
};

template <class TraitsT, class Visitor> 
std::ostream &operator<<(std::ostream &out, 
			 const Delaunay_triangulation_base_3<TraitsT, Visitor> &dt){
  dt.write(out);
  return out;
}

CGAL_KDS_END_INTERNAL_NAMESPACE

#endif

